
import IncrementalInference: selectFactorType, getDomain

export getDomain


## ============================================================================
# Starting integration with Manifolds.jl, via ApproxManifoldProducts.jl first
## ============================================================================

getDomain(::InstanceType{Point2Point2}) = Point2
getDomain(::InstanceType{Pose2Point2}) = Point2
getDomain(::InstanceType{Pose2Pose2}) = Pose2
# getDomain(::InstanceType{Pose3Point3}) = Point3
getDomain(::InstanceType{Pose3Pose3}) = Pose3
getDomain(::InstanceType{Pose2Point2BearingRange}) = BearingRange2

# defines on which variable manifold the measurements for a factor type are represented
getManifolds(fctType::InstanceType{Pose2Pose2}) = getManifolds(getDomain(fctType))
getManifolds(fctType::InstanceType{Pose2Point2}) = getManifolds(getDomain(fctType))
getManifolds(fctType::InstanceType{Point2Point2}) = getManifolds(getDomain(fctType))
getManifolds(fctType::InstanceType{Pose2Point2BearingRange}) = getManifolds(getDomain(fctType))



function selectFactorType(T1::Type{<:InferenceVariable}, T2::Type{<:InferenceVariable})
  # hacky version
  if T1 == Pose2 && T2 == Pose2
    return Pose2Pose2
  elseif T1 == Pose2 && T2 == Point2
    return Pose2Point2  
  elseif T1 == DynPose2 && T2 == DynPose2
      return DynPose2DynPose2
  elseif T1 == Point2 && T2 == Point2
    return Point2Point2
  elseif T1 == Point3 && T2 == Point3
    return Point3Point3
  elseif T1 == Pose3 && T2 == Pose3
    return Pose3Pose3
  else
    error("dont know which Factor type to select between $T1 and $T2")
  end
end




