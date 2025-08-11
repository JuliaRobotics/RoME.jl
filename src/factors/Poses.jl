
function (cf::CalcFactor{<:PriorPose2})(
    _m::AbstractArray{MT},
    _p::AbstractArray{PT},
) where {MT <: Real, PT <: Real}
    T = promote_type(MT, PT)
    m = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _m)
    p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
    return cf(m, p)
end

function (cf::CalcFactor{<:PriorPose2})(
    m::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
    p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
) where {T <: Real}
    M = getManifold(PriorPose2)
    X = log(M, p, m) # Currently X ∈ TₚM, #TODO should it be TₘM?
    return vee(M, p, X)
end

## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
compare(a::PriorPose2, b::PriorPose2; tol::Float64 = 1e-10) = compareDensity(a.Z, b.Z)

function (cf::CalcFactor{<:Pose2Pose2})(X, p, q)
  # X ∈ TₚM, X̂ ∈ TₚM, p,q ∈ M
  M = getManifold(Pose2)
  X̂ = log(M, p, q)
  return vee(M, p, X - X̂) # TODO check sign
end

# FIXME, rather have separate compareDensity functions
compare(a::Pose2Pose2, b::Pose2Pose2; tol::Float64 = 1e-10) = compareDensity(a.Z, b.Z)

# regular prior for Pose3

function (cf::CalcFactor{<:PriorPose3})(m, p)
    M = getManifold(Pose3)
    Xc = vee(M, p, log(M, p, m))
    return Xc
end

# Pose3Pose3 evaluation functions

function (cf::CalcFactor{<:Pose3Pose3})(X, p::ArrayPartition{T}, q) where {T}
    M = getManifold(Pose3Pose3)
    X̂ = log(M, p, q)
    Xc::SVector{6, T} = vee(M, p, X - X̂)
    return Xc
end

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
