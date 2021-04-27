
import IncrementalInference: selectFactorType
# import ApproxManifoldProducts: getManifold



## ============================================================================
# Starting integration with Manifolds.jl, via ApproxManifoldProducts.jl first
## ============================================================================

getManifold(::InstanceType{Point2Point2}) = Point2 |> getManifold
getManifold(::InstanceType{Pose2Point2}) = Point2 |> getManifold
getManifold(::InstanceType{Pose2Pose2}) = Pose2 |> getManifold
# getManifold(::InstanceType{Pose3Point3}) = Point3
getManifold(::InstanceType{Pose3Pose3}) = Pose3 |> getManifold
getManifold(::InstanceType{Pose2Point2BearingRange}) = BearingRange2 |> getManifold



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




