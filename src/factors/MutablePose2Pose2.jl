
export MutablePose2Pose2Gaussian, PackedMutablePose2Pose2Gaussian
# import IncrementalInference: getSample


"""
    $TYPEDEF

Specialized Pose2Pose2 factor type (Gaussian), which allows for rapid accumulation of odometry information as a branch on the factor graph.
"""
Base.@kwdef mutable struct MutablePose2Pose2Gaussian  <: IIF.AbstractManifoldMinimize
  Z::MvNormal
  timestamp::DateTime = now()
end
MutablePose2Pose2Gaussian(Z::MvNormal) = MutablePose2Pose2Gaussian(;Z)

DFG.getManifold(::MutablePose2Pose2Gaussian) = getManifold(Pose2) # Manifolds.SpecialEuclidean(2)


"""
    $SIGNATURES

Residual function for MutablePose2Pose2Gaussian.

Related

Pose2Pose2, Pose3Pose3, InertialPose3, DynPose2Pose2, Point2Point2, VelPoint2VelPoint2
"""
function (cf::CalcFactor{<:MutablePose2Pose2Gaussian})(X, p, q)
  M = getManifold(Pose2)
  q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
  #TODO allocalte for vee! see Manifolds #412, fix for AD
  Xc = zeros(3)
  vee!(M, Xc, q, log(M, q, q̂))
  return Xc
end


## Serialization support

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedMutablePose2Pose2Gaussian  <: AbstractPackedFactor
  Z::PackedSamplableBelief
  timestamp::Int64 # serialized in millisecond
end
function convert(::Type{MutablePose2Pose2Gaussian}, d::PackedMutablePose2Pose2Gaussian)
  return MutablePose2Pose2Gaussian(convert(SamplableBelief, d.datastr), unix2datetime(d.timestamp*1e-3))
end
function convert(::Type{PackedMutablePose2Pose2Gaussian}, d::MutablePose2Pose2Gaussian)
  return PackedMutablePose2Pose2Gaussian(convert(PackedSamplableBelief, d.Z), datetime2unix(d.timestamp)*1e3 |> Int)
end
