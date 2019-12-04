
export MutablePose2Pose2Gaussian, PackedMutablePose2Pose2Gaussian
# import IncrementalInference: getSample


"""
    $TYPEDEF

Specialized Pose2Pose2 factor type (Gaussian), which allows for rapid accumulation of odometry information as a branch on the factor graph.
"""
mutable struct MutablePose2Pose2Gaussian  <: IIF.FunctorPairwise
  Zij::MvNormal
end
function getSample(fct::MutablePose2Pose2Gaussian, N::Int=100)
  return (rand(fct.Zij, N), )
end

"""
    $SIGNATURES

Residual function for MutablePose2Pose2Gaussian.

Related

Pose2Pose2, Pose3Pose3, InertialPose3, DynPose2Pose2, Point2Point2, VelPoint2VelPoint2
"""
function (::MutablePose2Pose2Gaussian)(
                 res::Vector{Float64},
                 userdata,
                 idx::Int,
                 meas::Tuple,
                 wxi::Array{Float64,2},
                 wxj::Array{Float64,2}  )::Nothing
  #
  wXjhat = SE2(wxi[1:3,idx])*SE2(meas[1][1:3,idx])
  jXjhat = SE2(wxj[1:3,idx]) \ wXjhat
  se2vee!(res, jXjhat)
  nothing
end





## Serialization support

"""
$(TYPEDEF)
"""
mutable struct PackedMutablePose2Pose2Gaussian  <: IIF.PackedInferenceType
  datastr::String
  PackedMutablePose2Pose2Gaussian() = new()
  PackedMutablePose2Pose2Gaussian(x::String) = new(x)
end
function convert(::Type{MutablePose2Pose2Gaussian}, d::PackedMutablePose2Pose2Gaussian)
  return MutablePose2Pose2Gaussian(extractdistribution(d.datastr))
end
function convert(::Type{PackedMutablePose2Pose2Gaussian}, d::MutablePose2Pose2Gaussian)
  return PackedMutablePose2Pose2Gaussian(string(d.Zij))
end
