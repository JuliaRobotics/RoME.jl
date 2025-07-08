
#TODO deprecate, only use LieGroups
@inline function _vee(
    ::typeof(SpecialEuclidean(2; vectors = HybridTangentRepresentation())),
    X::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
) where {T <: Real}
    return SVector{3, T}(X.x[1][1], X.x[1][2], X.x[2][2])
end

#TODO deprecate, only use LieGroups
@inline function _compose(
    ::typeof(SpecialEuclidean(2; vectors = HybridTangentRepresentation())),
    p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
    q::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
) where {T <: Real}
    return ArrayPartition(p.x[1] + p.x[2] * q.x[1], p.x[2] * q.x[2])
end

function (cf::CalcFactor{<:PriorPose2})(
    _m::AbstractArray{MT},
    _p::AbstractArray{PT},
) where {MT <: Real, PT <: Real}
    T = promote_type(MT, PT)
    m = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _m)
    p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
    return cf(m, p)
end

# TODO the log here looks wrong (for gradients), consider:
# X = log(p⁻¹ ∘ m) 
# X = log(M, ϵ, Manifolds.compose(M, inv(M, p), m))
function (cf::CalcFactor{<:PriorPose2})(
    m::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
    p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
) where {T <: Real}
    M = getManifold(Pose2)
    ϵ = getPointIdentity(M)
    Xc = _vee(M, log(M, p, m))
    # X = log(M, ϵ, Manifolds.compose(M, inv(M, p), m))
    # Xc = vee(M, ϵ, X)
    return Xc
end

## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
compare(a::PriorPose2, b::PriorPose2; tol::Float64 = 1e-10) = compareDensity(a.Z, b.Z)

# Assumes X is a tangent vector
function (cf::CalcFactor{<:Pose2Pose2})(
    _X::AbstractArray{MT},
    _p::AbstractArray{PT},
    _q::AbstractArray{LT},
) where {MT, PT, LT}
    #TODO remove this convertions
    # @warn "This warning should not be triggered after StaticArrays upgrade" maxlog=10
    T = promote_type(MT, PT, LT)
    X = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _X)
    p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
    q = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _q)
    return cf(X, p, q)
end

# function calcPose2Pose2(
function (cf::CalcFactor{<:Pose2Pose2})(
    X::ArrayPartition{XT, Tuple{SVector{2, XT}, SMatrix{2, 2, XT, 4}}},
    p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
    q::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
) where {XT <: Real, T <: Real}
    M = getManifold(Pose2)
    # ϵ0 = ArrayPartition(zeros(SVector{2,T}), SMatrix{2, 2, T}(I))
    ϵ0 = getPointIdentity(M)

    ϵX = exp(M, ϵ0, X)
    # q̂ = Manifolds.compose(M, p, ϵX)    
    q̂ = _compose(M, p, ϵX)
    X_hat = log(M, q, q̂)#::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}
    # Xc = vee(M, q, X_hat)
    Xc = _vee(M, X_hat)#::SVector{3,T}
    return Xc
end

# FIXME, rather have separate compareDensity functions
compare(a::Pose2Pose2, b::Pose2Pose2; tol::Float64 = 1e-10) = compareDensity(a.Z, b.Z)

#

# regular prior for Pose3

function (cf::CalcFactor{<:PriorPose3})(m, p)
    M = getManifold(Pose3)
    Xc = vee(M, p, log(M, p, m))
    return Xc
end

# Pose3Pose3 evaluation functions

function (cf::CalcFactor{<:Pose3Pose3})(X, p::ArrayPartition{T}, q) where {T}
    M = getManifold(Pose3)
    # X: p_Xq̂
    # p: w_H_p
    # q: w_H_q
    # ̂q: w_H_̂q = w_H_p * exp_0(p_Xq̂)
    # w_H_̂q = Manifolds.compose(M, w_H_p, exp(M, getPointIdentity(M), p_Xq̂))

    q̂ = Manifolds.compose(M, p, exp(M, getPointIdentity(M), X))

    Xc::SVector{6, T} = get_coordinates(M, q, log(M, q, q̂), DefaultOrthogonalBasis())
    return Xc
end

# function (cf::CalcFactor{<:Pose3Pose3})(X, p, q)  
#   M = cf.manifold # getManifold(Pose3)
#   ϵX = exp(M, getPointIdentity(M), X)
#   q̂ = ArrayPartition(p.x[2]*ϵX.x[1] + p.x[1], p.x[2]*ϵX.x[2])
#   Xc = vee(M, q, log(M, q, q̂))
#   return Xc
# end

##
#TODO is this manifold not SO3
DFG.@defObservationType Pose3Pose3RotOffset AbstractManifoldMinimize Manifolds.SpecialEuclidean(3; vectors = HybridTangentRepresentation())

# measurement is in frame a, for example imu frame
# p and q is in frame b, for example body frame
# bRa is the rotation to get a in the b frame 
# measurement in frame a is converted to frame b and used to calculate the error
function (cf::CalcFactor{<:Pose3Pose3RotOffset})(aX, p, q, bRa)
    M = getManifold(Pose3Pose3RotOffset)
    # measurement in frame a, input is tangent, can also use vector transport 
    a_m = exp(M, getPointIdentity(M), aX)
    b_m = ArrayPartition(a_m.x[1], bRa * a_m.x[2])

    q̂ = Manifolds.compose(M, p, b_m)
    return vee(M, q, log(M, q, q̂)) # coordinates
end

##
DFG.@defObservationType Pose3Pose3Transform AbstractManifoldMinimize Manifolds.SpecialEuclidean(3; vectors = HybridTangentRepresentation())

function (cf::CalcFactor{<:Pose3Pose3Transform})(p_NX, p, q, Δ)
    M = getManifold(Pose3Pose3Transform)
    ε = getPointIdentity(M)

    Δn = compose(M, Δ, exp(M, ε, p_NX))
    q̂ = Manifolds.compose(M, p, Δn)

    Xc::SVector{6, T} = get_coordinates(M, q, log(M, q, q̂), DefaultOrthogonalBasis())
    return Xc
end

## ====================================
## Pose3Pose3UnitTrans Factor 

"""
  $(TYPEDEF)
Pose3Pose3 factor where the translation scale is not known, ie. Pose3Pose3 with unit (normalized) translation.
"""
DFG.@defObservationType Pose3Pose3UnitTrans AbstractManifoldMinimize Manifolds.SpecialEuclidean(3; vectors = HybridTangentRepresentation())

function (cf::CalcFactor{<:Pose3Pose3UnitTrans})(X, p::ArrayPartition{T}, q) where {T}
    M = getManifold(Pose3)
    q̂ = Manifolds.compose(M, p, exp(M, getPointIdentity(M), X))
    Xc::SVector{6, T} = get_coordinates(M, q, log(M, q, q̂), DefaultOrthogonalBasis())
    return SVector{6, T}(normalize(Xc[1:3])..., Xc[4:6]...)
end
