
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
  elseif T1 == Point2 && T2 == Point2
    return Point2Point2
  else
    error("dont know which Factor type to select between $T1 and $T2")
  end
end

buildFactorDefault(::Type{Pose2Pose2}) = Pose2Pose2(MvNormal(zeros(3), diagm([0.01;0.01;0.01])))
buildFactorDefault(::Type{Pose2Point2}) = Pose2Point2(MvNormal(zeros(2), diagm([0.01;0.01])))
buildFactorDefault(::Type{Pose2Point2BearingRange}) = Pose2Point2BearingRange(Normal(0,0.1),Normal(1,0.1))
buildFactorDefault(::Type{Point2Point2}) = Point2Point2(MvNormal(zeros(2), diagm([0.01;0.01])))
