using Manifolds
using StaticArrays
using Rotations
using LinearAlgebra
using DistributedFactorGraphs
using Dates

# FIXME use ApproxManifoldsProducts version instead
function TransformUtils.skew(v::SVector{3,T}) where T<:Real
    return SMatrix{3,3,T}(
            0,
         v[3],
        -v[2],
        -v[3],  
            0,  
         v[1],
         v[2],
        -v[1],
            0
    )
end

struct IMUDeltaManifold <: AbstractManifold{ℝ} end

# NOTE Manifold in not defined as a ProductManifold since we do not use the product metric. #701
# also see related SE₂(3) 

const IMUDeltaGroup = GroupManifold{ℝ, IMUDeltaManifold, MultiplicationOperation}

IMUDeltaGroup() = GroupManifold(IMUDeltaManifold(), MultiplicationOperation())

Manifolds.manifold_dimension(::IMUDeltaManifold) = 9

# Affine representation 
# Δ = [ΔR Δv Δp;
#      0   1 Δt;
#      0   0  1] ⊂ \R^5x5

#  ArrayPartition representation (TODO maybe swop order to [Δp; Δv; ΔR; Δt])
# Δ = [ΔR; Δv; Δp; Δt] 

function Manifolds.identity_element(M::IMUDeltaGroup) # was #SMatrix{5,5,Float64}(I)
    ArrayPartition(
        SMatrix{3,3,Float64}(I), # ΔR
        @SVector(zeros(3)),      # Δv
        @SVector(zeros(3)),      # Δp
        0.0                      # Δt
    )
end

DFG.getPointIdentity(M::IMUDeltaGroup) = identity_element(M)

function Manifolds.affine_matrix(G::IMUDeltaGroup, p::ArrayPartition{T}) where T<:Real
    return vcat(
        hcat(p.x[1], p.x[2], p.x[3]), 
        @SMatrix [0 0 0 1 p.x[4];
                  0 0 0 0 1]
    )
end

function vector_affine_matrix(G::IMUDeltaGroup, X::ArrayPartition{T}) where T<:Real
    return vcat(
        hcat(X.x[1], X.x[2], X.x[3]), 
        @SMatrix [0 0 0 0 X.x[4];
                  0 0 0 0 0]
    )
end

function Manifolds.inv(M::IMUDeltaGroup, p)
    ΔR = p.x[1]
    Δv = p.x[2]
    Δp = p.x[3]
    Δt = p.x[4]

    return ArrayPartition(
         ΔR',   
        -ΔR'*Δv,
        -ΔR'*(Δp - Δv*Δt),
        -Δt
    )
end

function Manifolds.compose(M::IMUDeltaGroup, p, q)
    ΔR = p.x[1]
    Δv = p.x[2]
    Δp = p.x[3]
    Δt = p.x[4]

    δR = q.x[1]
    δv = q.x[2]
    δp = q.x[3]
    δt = q.x[4]

    return ArrayPartition(
        ΔR * δR,
        Δv + ΔR * δv,
        Δp + Δv * δt + ΔR * δp,
        Δt + δt
    )
end

function Manifolds.hat(M::IMUDeltaGroup, Xⁱ::SVector{10, T}) where T<:Real
    return ArrayPartition(
        skew(Xⁱ[SA[7:9...]]), # θ ωΔt
        Xⁱ[SA[4:6...]],       # ν aΔt
        Xⁱ[SA[1:3...]],       # ρ vΔt
        Xⁱ[10],               # Δt
    )
end

function Manifolds.vee(M::IMUDeltaGroup, X::ArrayPartition{T}) where T<:Real
    return SVector{10,T}(
        X.x[3]...,   # ν aΔt 4:6
        X.x[2]...,   # ρ vΔt 1:3
        X.x[1][3,2], # θ⃗ₓ[3,2] 7
        X.x[1][1,3], # θ⃗ₓ[1,3] 8 
        X.x[1][2,1], # θ⃗ₓ[2,1] 9
        X.x[4], # Δt 10
    )
end

# small angle approx:
# Q ≈ I
# P ≈ 1/2 * I

function _Q(θ⃗)
    T = eltype(θ⃗)
    θ = norm(θ⃗)
    if θ ≈ 0
        return SMatrix{3,3,T}(I)
    else
        u = θ⃗/θ
        sθ, cθ = sincos(θ)
        uₓ = skew(u)
        # NOTE difference in references here --- (θ - sθ)/θ^2 vs (θ - sθ)/θ
        # with no ^2 looking correct when compared to exp of SE3
        return SMatrix{3,3,T}(I) + (1 - cθ)/θ * uₓ + (θ - sθ)/θ * uₓ^2
    end
end

function _P(θ⃗)
    T = eltype(θ⃗)
    θ = norm(θ⃗)
    if θ ≈ 0
        return 1/2*SMatrix{3,3,T}(I)
    else
        u = θ⃗/θ
        sθ, cθ = sincos(θ)
        uₓ = skew(u)
        return 1/2*SMatrix{3,3,T}(I) + (θ - sθ)/θ^2 * uₓ + (cθ + 1/2*θ^2 - 1)/θ^2 * uₓ^2
    end
end

#TODO rename to exp_lie?
function Manifolds.exp(M::IMUDeltaGroup, X::ArrayPartition{T}) where T<:Real
    θ⃗ₓ = X.x[1] # ωΔt
    
    ν = X.x[2]  # aΔt
    ρ = X.x[3]  # vΔt

    Δt = X.x[4]

    # ωΔt = vee(θ⃗ₓ)
    θ⃗ = SA[θ⃗ₓ[3,2]; θ⃗ₓ[1,3]; θ⃗ₓ[2,1]]

    P = _P(θ⃗)
    Q = _Q(θ⃗)

    M_SO3 = SpecialOrthogonal(3)
    q = ArrayPartition(
        exp(M_SO3, getPointIdentity(M_SO3), θ⃗ₓ),
        Q*ν,
        Q*ρ + P*ν*Δt,
        Δt,
    )
    return q
end

function Manifolds.exp(M::IMUDeltaGroup, p::ArrayPartition{T}, X::ArrayPartition{T}) where T<:Real
    q = exp(M, X)
    return Manifolds.compose(M, p, q)
end

function Manifolds.log(M::IMUDeltaGroup, p)
    ΔR = p.x[1]
    Δv = p.x[2]
    Δp = p.x[3]
    Δt = p.x[4]

    Rv = RotationVec(ΔR)
    θ⃗ = SA[Rv.sx, Rv.sy, Rv.sz]

    P = _P(θ⃗)
    Q = _Q(θ⃗)
    iQ = inv(Q)
    
    return ArrayPartition(
        log_lie(SpecialOrthogonal(3), ΔR), # θ⃗ₓ
        iQ*Δv, # ν aΔt
        iQ*(Δp - P*iQ*Δv*Δt), # ρ vΔt 
        Δt
    )
end

#TODO TEST 
function Manifolds.log(M::IMUDeltaGroup, p, q)
    return log(M, Manifolds.compose(M, inv(M, p), q))
end


# TODO consolidate with Manifolds notation 
# function Manifolds.log(M::IMUDeltaGroup, p::SMatrix{5,5,T}, q::SMatrix{5,5,T})
#     qinvp = Manifolds.compose(M, inv(M, q), p)
#     .....
# end

# compute the expected delta from p to q on the IMUDeltaGroup
# ⊟ = boxminus
function boxminus(::IMUDeltaGroup, p, q; g⃗ = SA[0,0,9.81]) 
    Rᵢ = p.x[1]
    vᵢ = p.x[2]
    pᵢ = p.x[3]
    tᵢ = p.x[4]

    Rⱼ = q.x[1]
    vⱼ = q.x[2]
    pⱼ = q.x[3]
    tⱼ = q.x[4]

    Δt = tⱼ - tᵢ
    ΔR = Rᵢ'Rⱼ

    Δp = Rᵢ' * (pⱼ - pᵢ - vᵢ * Δt + 1/2 * g⃗ * Δt^2)
    Δv = Rᵢ' * (vⱼ - vᵢ + g⃗ * Δt)

    return ArrayPartition(
        ΔR,
        Δv,
        Δp,
        Δt
    )
end

#TODO test, unused, likeley to be replaced because of ambiguity
# right-⊖ : Xₚ = q ⊖ p = log(p⁻¹∘q)
function ominus(M::IMUDeltaGroup, q, p)
    log(M, Manifolds.compose(M, inv(M, p), q))
end

# small adjoint ad
function adjointMatrix(::IMUDeltaGroup, X::ArrayPartition{T}) where T
    θ⃗ₓ = X.x[1] # ωΔt
    ν = X.x[2]  # aΔt
    ρ = X.x[3]  # vΔt

    IΔt = SMatrix{3,3,T}(X.x[4]*I)

    ρₓ = skew(ρ)
    νₓ = skew(ν)
    
    m0 = zeros(3,3) 
    v0 = zeros(3)

    return SMatrix{10,10,T}([
         θ⃗ₓ  -IΔt  ρₓ   ν;
         m0    θ⃗ₓ  νₓ  v0;
         m0    m0  θ⃗ₓ  v0;
              zeros(1,10);
    ])

end

# Adjoint Ad
function AdjointMatrix(::IMUDeltaGroup, p::ArrayPartition{T}) where T
    ΔR = p.x[1]
    Δv = p.x[2]
    Δp = p.x[3]
    Δt = p.x[4]

    Δvₓ = skew(Δv)
    pmvtₓ = skew(Δp - Δv*Δt)

    m0 = zeros(3,3) 
    v0 = zeros(3)

    return SMatrix{10,10,T}([
        ΔR  -ΔR*Δt  pmvtₓ*ΔR  Δv;
        m0      ΔR    Δvₓ*ΔR  v0;
        m0      m0        ΔR  v0;
        zeros(1,9)             1;
    ])

end

# right Jacobian
# FIXME moving general Lie Group Jacobian up to ApproxManifoldProducts
function Jr(M::IMUDeltaGroup, X; order=5)
    adx = adjointMatrix(M, X)
    mapreduce(+, 0:order) do i
        (-adx)^i / factorial(i + 1)
    end
end

## ======================================================================================
## IMU DELTA FACTOR 
## ======================================================================================
struct IMUMeasurements
    accelerometer::Vector{SVector{3,Float64}}
    gyroscope::Vector{SVector{3,Float64}}
    timestamps::Vector{Float64}
    Σy::SMatrix{6,6,Float64}
end

"""
$TYPEDEF

Factor type for inertial odometry (preintegration).
"""
Base.@kwdef struct IMUDeltaFactor{T <: SamplableBelief} <: AbstractManifoldMinimize
    Z::T # NOTE dim is 9 as Δt is not included in covariance
    Δt::Float64
    Δ::ArrayPartition{Float64, Tuple{SMatrix{3, 3, Float64, 9}, SVector{3, Float64}, SVector{3, Float64}, Float64}}
    Σ::SMatrix{10,10,Float64}
    J_b::SMatrix{10,6,Float64} = zeros(SMatrix{10,6,Float64})
    # accelerometer bias, gyroscope bias 
    b::SVector{6, Float64} = zeros(SVector{6, Float64})
    #optional raw measurements
    raw_measurements::Union{Nothing,IMUMeasurements} = nothing
end

function IIF.getSample(cf::CalcFactor{<:IMUDeltaFactor})
    return exp(IMUDeltaGroup(), hat(IMUDeltaGroup(), SA[rand(cf.factor.Z)..., cf.factor.Δt]))
end

function IIF.getFactorMeasurementParametric(f::IMUDeltaFactor)
    iΣ = convert(SMatrix{9,9,Float64,81}, invcov(f.Z))
    return f.Δ, iΣ
end

IIF.getManifold(::IMUDeltaFactor) = IMUDeltaGroup()

function IIF.preambleCache(fg::AbstractDFG, vars::AbstractVector{<:DFGVariable}, ::IMUDeltaFactor)
    (timestams=(vars[1].nstime,vars[2].nstime),)
end

# factor residual

function (cf::CalcFactor{<:IMUDeltaFactor})(
    Δmeas, # point on IMUDeltaGroup
    p::ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}, SVector{3, T}, T}}, 
    q::ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}, SVector{3, T}, T}},
    b::SVector{6,T} = zeros(SVector{6,T})
) where T <: Real
    #
    M = IMUDeltaGroup()
    # imu measurment Delta, corrected for bias with # b̄ = cf.factor.b
    Δi = Manifolds.compose(M, Δmeas, exp(M, hat(M, cf.factor.J_b * (b - cf.factor.b))))
    # expected delta from p to q
    Δhat = boxminus(M, p, q)
    # residual 
    Xhat = log(M, Manifolds.compose(M, inv(M, Δi), Δhat))
    # should not include Δt only 1:9
    return vee(M, Xhat)[SOneTo(9)]
end

function (cf::CalcFactor{<:IMUDeltaFactor})(
    Δmeas, # point on IMUDeltaGroup
    _p::ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}, SVector{3, T}}}, 
    _q::ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}, SVector{3, T}}},
    b::SVector{6,T} = zeros(SVector{6,Float64})
) where T <: Real
    p_t = Dates.value(cf.cache.timestams[1])*1e-9
    q_t = Dates.value(cf.cache.timestams[2])*1e-9
    p = ArrayPartition(_p.x[1], _p.x[2], _p.x[3], p_t)
    q = ArrayPartition(_q.x[1], _q.x[2], _q.x[3], q_t)
    return cf(Δmeas, p, q, b)
end

# ArrayPartition{Float64, Tuple{MMatrix{3, 3, Float64, 9}, MVector{3, Float64}, MVector{3, Float64}}}
function (cf::CalcFactor{<:IMUDeltaFactor})(
    Δmeas, # point on IMUDeltaGroup
    _p::ArrayPartition{Float64, Tuple{MMatrix{3, 3, Float64, 9}, MVector{3, Float64}, MVector{3, Float64}}}, 
    _q::ArrayPartition{Float64, Tuple{MMatrix{3, 3, Float64, 9}, MVector{3, Float64}, MVector{3, Float64}}},
    b::AbstractVector = zeros(SVector{6,Float64})
)
    p_t = Dates.value(cf.cache.timestams[1])*1e-9
    q_t = Dates.value(cf.cache.timestams[2])*1e-9
    p = ArrayPartition(SMatrix{3,3,Float64,9}(_p.x[1]), SVector{3}(_p.x[2]), SVector{3}(_p.x[3]), p_t)
    q = ArrayPartition(SMatrix{3,3,Float64,9}(_q.x[1]), SVector{3}(_q.x[2]), SVector{3}(_q.x[3]), q_t)
    return cf(Δmeas, p, q, SVector{6}(b))
end

function (cf::CalcFactor{<:IMUDeltaFactor})(
    Δmeas, 
    p_SE3,
    p_vel, 
    q_SE3,
    q_vel,
    b::SVector{6,T} = zeros(SVector{6,Float64})
) where T <: Real
    p_t = Dates.value(cf.cache.timestams[1])*1e-9
    q_t = Dates.value(cf.cache.timestams[2])*1e-9
    p = ArrayPartition(p_SE3.x[2], p_vel, p_SE3.x[1], p_t)
    q = ArrayPartition(q_SE3.x[2], q_vel, q_SE3.x[1], q_t)
    return cf(Δmeas, p, q, b)
end

function _τδt(δt)
    δtI = SMatrix{3,3,Float64}(δt*I)
    m = zeros(10,6)
    m[4:6,1:3] .= δtI
    m[7:9,4:6] .= δtI
    return SMatrix{10,6,Float64}(m)
end

function integrateIMUDelta(Δij, Σij, Δij_J_b, a, ω, a_b, ω_b, δt, Σy)
    M = IMUDeltaGroup()
    ωδt = (ω - ω_b) * δt
    aδt = (a - a_b) * δt
    Xc = SVector{10,Float64}(vcat([0., 0, 0], aδt, ωδt, δt))

    X = hat(M, Xc)
    δjk = exp(M, X)

    Δik = Manifolds.compose(M, Δij, δjk)

    # Jacobians
    τ_J_y = _τδt(δt)
    τ_J_b = -_τδt(δt)
    δ_J_τ = Jr(M, X)

    Δik_J_Δij = inv(AdjointMatrix(M, δjk))
    #? Δik_J_Δij = AdjointMatrix(M, inv(M, δjk))
    Δik_J_δjk = I # 10x10

    # Propagate Covariance
    Δik_J_y = Δik_J_δjk * δ_J_τ * τ_J_y
    Σik = Δik_J_Δij * Σij * Δik_J_Δij' + Δik_J_y * Σy * Δik_J_y'

    # Jacobian wrt biases
    Δik_J_b = Δik_J_Δij * Δij_J_b  + Δik_J_δjk * δ_J_τ * τ_J_b

    return Δik, Σik, Δik_J_b
end

#assuming accels and gyros same rate here
function preintegrateIMU(accels, gyros, timestamps, Σy, a_b, ω_b)
    M = IMUDeltaGroup()
    Δ   = identity_element(M)
    Σ   = zeros(10,10)
    J_b = zeros(10,6)
    
    for i in eachindex(timestamps)[2:end]
        a = accels[i]
        ω = gyros[i]    
        dt = timestamps[i] - timestamps[i-1]
        Δ, Σ, J_b = integrateIMUDelta(Δ, Σ, J_b, a, ω, a_b, ω_b, dt, Σy)
    end
    Δ, Σ, J_b
end

function IMUDeltaFactor(
    accels::AbstractVector,
    gyros::AbstractVector,
    timestamps::AbstractVector,
    Σy,
    a_b = SA[0.,0,0],
    ω_b = SA[0.,0,0],
)
    M = IMUDeltaGroup()
    Δ, Σ, J_b = preintegrateIMU(accels, gyros, timestamps, Σy, a_b, ω_b)
    Δt = Δ.x[4]
    
    Xc = vee(M, log(M, Δ))
    
    SM = SymmetricPositiveDefinite(9)
    S = Σ[1:9,1:9]
    ch = check_point(SM, S; atol = 1e-9) 
    if !isnothing(ch)
        @warn "IMU Covar check" ch
        S = (S + S') / 2
        S = S + diagm((diag(S) .== 0)*1e-15)
        ch = check_point(SM, S)
        !isnothing(ch) && @error "IMU Covar check" ch
    end

    Z = MvNormal(Xc[1:9], S)

    IMUDeltaFactor(
        Z,
        Δt,
        Δ,
        SMatrix{10,10,Float64}(Σ),
        J_b,
        SA[a_b...; ω_b...],
        IMUMeasurements(
            accels,
            gyros,
            timestamps,
            Σy
        )
    )
end
