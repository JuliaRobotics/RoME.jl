using Manifolds
using StaticArrays
using Rotations
using LinearAlgebra
using DistributedFactorGraphs
using Dates


struct SpecialGalileanManifold <: AbstractManifold{ℝ} end

# NOTE Manifold in not defined as a ProductManifold since we do not use the product metric. #701
# also see related SE₂(3) 
# SemidirectProductGroup(SpecialEuclidean(N), (ProductGroup(TranslationGroup(N), TranslationGroup(1))))
# (SO(n) ⋉ ℝⁿ) ⋉ (ℝⁿ × ℝ)

"""
    SpecialGalileanGroup

References: 
- https://hal.science/hal-02183498/document
- TODO new reference: https://arxiv.org/pdf/2312.07555
- TODO new reference: https://arxiv.org/pdf/2409.14276

Affine representation 
Δ = [ΔR Δv Δp;
     0   1 Δt;
     0   0  1] ⊂ ℝ⁵ˣ⁵

ArrayPartition representation (TODO maybe swop order to [Δp; Δv; ΔR; Δt])
Δ = [ΔR; Δv; Δp; Δt] 
"""
const SpecialGalileanGroup = GroupManifold{ℝ, SpecialGalileanManifold, MultiplicationOperation}

SpecialGalileanGroup() = GroupManifold(SpecialGalileanManifold(), MultiplicationOperation(), LeftInvariantRepresentation())

Manifolds.manifold_dimension(::SpecialGalileanManifold) = 9


function Manifolds.identity_element(M::SpecialGalileanGroup) # was #SMatrix{5,5,Float64}(I)
    ArrayPartition(
        SMatrix{3,3,Float64}(I), # ΔR
        @SVector(zeros(3)),      # Δv
        @SVector(zeros(3)),      # Δp
        0.0                      # Δt
    )
end

DFG.getPointIdentity(M::SpecialGalileanGroup) = identity_element(M)

function Manifolds.affine_matrix(G::SpecialGalileanGroup, p::ArrayPartition{T}) where T<:Real
    return vcat(
        hcat(p.x[1], p.x[2], p.x[3]), 
        @SMatrix [0 0 0 1 p.x[4];
                  0 0 0 0 1]
    )
end

function vector_affine_matrix(G::SpecialGalileanGroup, X::ArrayPartition{T}) where T<:Real
    return vcat(
        hcat(X.x[1], X.x[2], X.x[3]), 
        @SMatrix [0 0 0 0 X.x[4];
                  0 0 0 0 0]
    )
end

function Manifolds.inv(M::SpecialGalileanGroup, p)
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

function Manifolds.compose(M::SpecialGalileanGroup, p, q)
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

function Manifolds.hat(M::SpecialGalileanGroup, Xⁱ::SVector{10, T}) where T<:Real
    return ArrayPartition(
        ApproxManifoldProducts.skew(Xⁱ[SA[7:9...]]), # θ ωΔt
        Xⁱ[SA[4:6...]],       # ν aΔt
        Xⁱ[SA[1:3...]],       # ρ vΔt
        Xⁱ[10],               # Δt
    )
end

function Manifolds.vee(M::SpecialGalileanGroup, X::ArrayPartition{T}) where T<:Real
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
        uₓ = ApproxManifoldProducts.skew(u)
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
        uₓ = ApproxManifoldProducts.skew(u)
        return 1/2*SMatrix{3,3,T}(I) + (θ - sθ)/θ^2 * uₓ + (cθ + 1/2*θ^2 - 1)/θ^2 * uₓ^2
    end
end

#TODO rename to exp_lie?
Manifolds.exp(M::SpecialGalileanGroup, X::ArrayPartition{T}) where T<:Real = error("use exp_lie instead")
function Manifolds.exp_lie(M::SpecialGalileanGroup, X::ArrayPartition{T}) where T<:Real
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

#TODO is this now exp_inv? to fit with Manifold.jl (until LieGroups.jl is done)
function Manifolds.exp(M::SpecialGalileanGroup, p::ArrayPartition{T}, X::ArrayPartition{T}) where T<:Real
    q = exp_lie(M, X)
    return Manifolds.compose(M, p, q)
end

Manifolds.log(M::SpecialGalileanGroup, p::ArrayPartition{T}) where T<:Real = error("use log_lie instead")
function Manifolds.log_lie(M::SpecialGalileanGroup, p)
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
#TODO is this now log_inv? to fit with Manifold.jl (until LieGroups.jl is done)
# right-⊖ : Xₚ = q ⊖ p = log(p⁻¹∘q)
function Manifolds.log(M::SpecialGalileanGroup, p, q)
    return log_lie(M, Manifolds.compose(M, inv(M, p), q))
end

# compute the expected delta from p to q on the SpecialGalileanGroup
# ⊟ = boxminus
function boxminus(::SpecialGalileanGroup, p, q; g⃗ = SA[0,0,9.81]) 
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

# small adjoint ad
function adjointMatrix(::SpecialGalileanGroup, X::ArrayPartition{T}) where T
    θ⃗ₓ = X.x[1] # ωΔt
    ν = X.x[2]  # aΔt
    ρ = X.x[3]  # vΔt

    IΔt = SMatrix{3,3,T}(X.x[4]*I)

    ρₓ = ApproxManifoldProducts.skew(ρ)
    νₓ = ApproxManifoldProducts.skew(ν)
    
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
function AdjointMatrix(::SpecialGalileanGroup, p::ArrayPartition{T}) where T
    ΔR = p.x[1]
    Δv = p.x[2]
    Δp = p.x[3]
    Δt = p.x[4]

    Δvₓ = ApproxManifoldProducts.skew(Δv)
    pmvtₓ = ApproxManifoldProducts.skew(Δp - Δv*Δt)

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
function Jr(M::SpecialGalileanGroup, X; order=5)
    adx = adjointMatrix(M, X)
    mapreduce(+, 0:order) do i
        (-adx)^i / factorial(i + 1)
    end
end

## ======================================================================================
## IMU DELTA FACTOR 
## ======================================================================================
struct IMUMeasurement
    accelerometer::SVector{3,Float64}
    gyroscope::SVector{3,Float64}
    deltatime::Float64
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
    #optional raw measurements, NOTE used for development and may be removed in the future
    raw_measurements::Vector{IMUMeasurement} = IMUMeasurement[]
end

function IIF.getSample(cf::CalcFactor{<:IMUDeltaFactor})
    return exp_lie(SpecialGalileanGroup(), hat(SpecialGalileanGroup(), SA[rand(cf.factor.Z)..., cf.factor.Δt]))
end

function IIF.getFactorMeasurementParametric(f::IMUDeltaFactor)
    iΣ = convert(SMatrix{9,9,Float64,81}, invcov(f.Z))
    return f.Δ, iΣ
end

IIF.getManifold(::IMUDeltaFactor) = SpecialGalileanGroup()

function IIF.preambleCache(fg::AbstractDFG, vars::AbstractVector{<:DFGVariable}, ::IMUDeltaFactor)
    if vars[1] isa DFGVariable{<:Pose3}
        (timestams = (vars[1].nstime, vars[3].nstime),)
    else
        (timestams = (vars[1].nstime, vars[2].nstime),)
    end

end

# factor residual

function (cf::CalcFactor{<:IMUDeltaFactor})(
    Δmeas, # point on SpecialGalileanGroup
    p::ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}, SVector{3, T}, T}}, 
    q::ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}, SVector{3, T}, T}},
    b::SVector{6,T} = zeros(SVector{6,T})
) where T <: Real
    #
    M = SpecialGalileanGroup()
    # imu measurment Delta, corrected for bias with # b̄ = cf.factor.b #TODO check if (b - cf.factor.b) is correct
    Δi = Manifolds.compose(M, Δmeas, exp_lie(M, hat(M, cf.factor.J_b * (b - cf.factor.b))))
    # expected delta from p to q
    Δhat = boxminus(M, p, q)
    # residual 
    Xhat = log_lie(M, Manifolds.compose(M, inv(M, Δi), Δhat))

    Xc_hat = vee(M, Xhat)
    @assert isapprox(Δi.x[4], Δhat.x[4], atol=1e-6) "Time descrepancy in IMUDeltaFactor: Δt = $(Xc_hat[10])"
    # should not include Δt only 1:9
    return Xc_hat[SOneTo(9)]
end

function (cf::CalcFactor{<:IMUDeltaFactor})(
    Δmeas, # point on SpecialGalileanGroup
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
    Δmeas, # point on SpecialGalileanGroup
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
    b = zeros(SVector{6,Float64})
) where T <: Real
    p = ArrayPartition(p_SE3.x[2], p_vel, p_SE3.x[1])
    q = ArrayPartition(q_SE3.x[2], q_vel, q_SE3.x[1])
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
    M = SpecialGalileanGroup()
    ωδt = (ω - ω_b) * δt
    aδt = (a - a_b) * δt
    Xc = SVector{10,Float64}(vcat([0., 0, 0], aδt, ωδt, δt))

    X = hat(M, Xc)
    δjk = exp_lie(M, X)

    Δik = Manifolds.compose(M, Δij, δjk)

    # Jacobians
    τ_J_y = _τδt(δt)
    τ_J_b = -_τδt(δt)
    δ_J_τ = Jr(M, X) #  jacobian(δjk=exp(X) wrt X=τ) = right jacobian

    # Σik = Δik_J_Δij * Σij * Δik_J_Δij'
    # this is jacobian(Δik = Δij ∘ δjk) wrt Δij = inv(Ad(δjk))
    # we can maybe test with
    # differentiate Δij∘δjk with respect to Δij in the direction X (tangent at Δij)
    # translate_diff(G, δjk, Δij, X, Manifolds.RightBackwardAction())
    # X will be the basis vectors to build up the jacobian
    Δik_J_Δij = inv(AdjointMatrix(M, δjk))
    #? Δik_J_Δij = AdjointMatrix(M, inv(M, δjk))
    Δik_J_δjk = I # 10x10 jacobian(Δik = Δij ∘ δjk) wrt δjk = I

    # Propagate Covariance
    Δik_J_y = Δik_J_δjk * δ_J_τ * τ_J_y
    Σik = Δik_J_Δij * Σij * Δik_J_Δij'  +  Δik_J_y * Σy * Δik_J_y'

    # Jacobian wrt biases
    Δik_J_b = Δik_J_Δij * Δij_J_b  +  Δik_J_δjk * δ_J_τ * τ_J_b

    return Δik, Σik, Δik_J_b
end

#assuming accels and gyros same rate here
function preintegrateIMU(accels, gyros, deltatimes, Σy, a_b, ω_b)
    M = SpecialGalileanGroup()
    Δ   = identity_element(M)
    Σ   = zeros(10,10)
    J_b = zeros(10,6)
    
    for (a, ω, dt) in zip(accels, gyros, deltatimes)
        Δ, Σ, J_b = integrateIMUDelta(Δ, Σ, J_b, a, ω, a_b, ω_b, dt, Σy)
    end
    Δ, Σ, J_b
end

function IMUDeltaFactor(
    accels::AbstractVector,
    gyros::AbstractVector,
    deltatimes::AbstractVector,
    Σy,
    a_b = SA[0.,0,0],
    ω_b = SA[0.,0,0],
)
    M = SpecialGalileanGroup()
    Δ, Σ, J_b = preintegrateIMU(accels, gyros, deltatimes, Σy, a_b, ω_b)
    Δt = Δ.x[4]
    
    Xc = vee(M, log_lie(M, Δ))
    
    SM = SymmetricPositiveDefinite(9)
    S = Σ[1:9,1:9]
    ch = check_point(SM, S; atol = 1e-9) 
    if !isnothing(ch)
        @warn "IMU Covar check" ch maxlog=1
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
        SMatrix{10,6,Float64}(J_b),
        SA[a_b...; ω_b...],
        IMUMeasurement[]
    )
end



## serde

struct PackedIMUDeltaFactor{T <: PackedSamplableBelief} <: AbstractPackedFactor
    Z::T # NOTE dim is 9 as Δt is not included in covariance
    dt::Float64
    D::Vector{Float64}
    Sigma::Vector{Float64} #SMatrix{10,10,Float64}
    J_b::Vector{Float64} # SMatrix{10,6,Float64}
    # accelerometer bias, gyroscope bias 
    b::Vector{Float64}
end

function PackedIMUDeltaFactor(;
    Z,
    dt,
    D,
    Sigma,
    J_b,
    b
)
    _gettype(zt::PackedSamplableBelief) = zt
    _gettype(zt) = DistributedFactorGraphs.getTypeFromSerializationModule(zt["_type"])(; zt...)

    _Z = _gettype(Z)
    # _Z = ZT(; Z...)
    _dt = Float64(dt)
    _D = Float64.(D)
    _Sigma = Float64.(Sigma)
    _J_b = Float64.(J_b)
    _b = Float64.(b)
    PackedIMUDeltaFactor(
        _Z,
        _dt,
        _D,
        _Sigma,
        _J_b,
        _b
    )
end

function convert(::Type{<:PackedIMUDeltaFactor}, d::IMUDeltaFactor)
    Z = convert(PackedSamplableBelief, d.Z)
    return PackedIMUDeltaFactor(;
        Z,
        dt = d.Δt,
        D = d.Δ[1:end],
        Sigma = collect(d.Σ[:]),
        J_b = collect(d.J_b[:]),
        b = collect(d.b),
    )
end

function convert(::Type{<:IMUDeltaFactor}, d::PackedIMUDeltaFactor)
    @assert 16 == length(d.D) "Deserializing a PackedIMUDeltaFactor not 16 in length, .D = $(length(d.D))"
    return IMUDeltaFactor(
        convert(SamplableBelief, d.Z),
        Float64(d.dt),
        ArrayPartition(
            SMatrix{3, 3, Float64, 9}(d.D[1:9]),
            SVector{3, Float64}(d.D[10:12]),
            SVector{3, Float64}(d.D[13:15]),
            Float64(d.D[16])
        ),
        SMatrix{10, 10, Float64}(d.Sigma),
        SMatrix{10,6,Float64}(d.J_b),
        SVector{6, Float64}(d.b),
        IMUMeasurement[],
    )
end