
import IncrementalInference: selectFactorType

## =============================================================================
## Needs a home
## =============================================================================

getManifolds(vartype::Type{Pose2Pose2}) = getManifolds(Pose2)
getManifolds(vartype::Type{Pose2Point2}) = getManifolds(Point2)
getManifolds(vartype::Type{Point2Point2}) = getManifolds(Point2)
getManifolds(vartype::Type{Pose2Point2BearingRange}) = (:Circular, :Euclid)


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
