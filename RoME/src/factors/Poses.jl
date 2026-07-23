# NOTE we follow `residual = measurement - prediction` as a convention

# NOTE ON PRIORS ON MANIFOLDS:
# For prior factors, we compute the residual as log(M, p, m), which gives a tangent vector at the current state estimate p
# pointing toward the measurement m. This means the residual lives in the tangent space at p (TₚM).
#
# - The optimizer linearizes and updates in the tangent space at the current estimate p.

# Note on injectivity_radius
# If the measurement is outside the injectivity radius of the manifold, the log will not be valid.
# even if we use Xq = log(M, q, q̂), we need to transport Xq to Xp, parallel_transport_to still uses the log map.

function (cf::CalcFactor{<:PriorPose2})(m, p)
    M = getManifold(PriorPose2)
    X = log(M, p, m) # the residual is calculated at the current state estimate p, ie. X ∈ TₚM
    return vee(LieAlgebra(M), X)
end

function (cf::CalcFactor{<:Pose2Pose2})(X, p, q)
    # X ∈ TₚM, X̂ ∈ TₚM, p,q ∈ M
    M = getManifold(Pose2Pose2)
    X̂ = log(M, p, q)
    return vee(LieAlgebra(M), X - X̂)
end

function (cf::CalcFactor{<:PriorPose3})(m, p)
    M = getManifold(PriorPose3)
    return vee(LieAlgebra(M), log(M, p, m))
end

function (cf::CalcFactor{<:Pose3Pose3})(X, p, q)
    # X ∈ TₚM, X̂ ∈ TₚM, p,q ∈ M
    M = getManifold(Pose3Pose3)
    X̂ = log(M, p, q)
    return vee(LieAlgebra(M), X - X̂)
end

# FIXME, rather have separate compareDensity functions
compare(a::Pose2Pose2, b::Pose2Pose2; tol::Float64 = 1e-10) = compareDensity(a.Z, b.Z)
## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
compare(a::PriorPose2, b::PriorPose2; tol::Float64 = 1e-10) = compareDensity(a.Z, b.Z)

##
DFG.@defObservationType Pose3Pose3RotOffset RelativeObservation LeftInvariantMetricSE(
    3,
)

# Measurement `Xs` is a tangent vector in frame 's' (sensor frame)
# `pRs` is the rotation offset to get frame 'p' into frame 's'
function (cf::CalcFactor{<:Pose3Pose3RotOffset})(Xs, p, q, pRs)
    # X̂p ∈ TₚM; p, q ∈ M
    M = getManifold(Pose3Pose3RotOffset)
    
    # X̂p ∈ TₚM is the relative tangent vector between the body poses, anchored at base point `p`.
    X̂p = log(M, p, q)
        
    # By evaluating `diff_left_compose` at base point `p`, we push the body's tangent vector `X̂p` 
    # forward to the point on the manifold where the sensor exists (p ∘ pTs).
    # This yields `X̂s` (anchored at p ∘ pTs), aligning with measurement `Xs`.
    pTs = ArrayPartition(zeros(eltype(pRs), 3), pRs)
    X̂s = diff_left_compose(base_lie_group(M), p, pTs, X̂p)
    
    # NOTE we could have used the adjoint as well to transform the predicted tangent vector to the SENSOR frame 's'
    # sRp = transpose(pRs)
    # sTp = ArrayPartition(zeros(eltype(pRs), 3), sRp)
    # X̂s = adjoint(base_lie_group(M), sTp, X̂p)

    # Calculate the residual
    return vee(LieAlgebra(M), Xs - X̂s)
end

#
DFG.@defObservationType Pose3Pose3Offset RelativeObservation LeftInvariantMetricSE(3)

function (cf::CalcFactor{<:Pose3Pose3Offset})(X, p, q, pTs)
    M = getManifold(Pose3Pose3Offset)
    X̂p = log(M, p, q)
    X̂ = diff_left_compose(base_lie_group(M), p, pTs, X̂p)
    return vee(LieAlgebra(M), X - X̂)
end

## ====================================
## Pose3Pose3UnitTrans Factor 

"""
  $(TYPEDEF)
Pose3Pose3 factor where the translation scale is not known, ie. Pose3Pose3 with unit (normalized) translation.
"""
DFG.@defObservationType Pose3Pose3UnitTrans RelativeObservation LeftInvariantMetricSE(3)

# NOTE: this manifold is actually Sphere(2) × SpecialOrthogonal(3), but we embed in LeftInvariantMetricSE(3) and use the chordal distance `-`.
function (cf::CalcFactor{<:Pose3Pose3UnitTrans})(X, p, q)
    M = getManifold(Pose3Pose3UnitTrans)
    
    X̂ = log(M, p, q)
    X_coords = vee(LieAlgebra(M), X)
    X̂_coords = vee(LieAlgebra(M), X̂)
    
    # Calculate the residual at p: measurement - prediction
    # We normalize the prediction's translation direction to match the unit measurement X
    return SVector{6}(
        (normalize(X_coords[1:3]) - normalize(X̂_coords[1:3]))...,
        (X_coords[4:6] - X̂_coords[4:6])...
    )
end

# #  FIXME needed until AMP#41 is done hopefully can be removed soon 🐛💥
# # Base.convert(::Type{<:Tuple}, ::typeof(LeftInvariantMetricSE(2))) = (:Euclid,:Euclid,:Circular)
# AMP._manifoldtuple(::typeof(LeftInvariantMetricSE(2))) = (:Euclid, :Euclid, :Circular)
# function AMP._manifoldtuple(::typeof(LeftInvariantMetricSE(3)))
#     return (:Euclid, :Euclid, :Euclid, :Circular, :Circular, :Circular)
# end
