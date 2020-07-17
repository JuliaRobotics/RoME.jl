# Graph analysis tools and support calculation tools

import IncrementalInference: selectFactorType, buildFactorDefault


function rangeErrMaxPoint2(fgl1::AbstractDFG, id1, fgl2::AbstractDFG ,id2)
  mv1 = getKDEMax(getVertKDE(fgl1,id1))
  mv2 = getKDEMax(getVertKDE(fgl2,id2))
  return norm(mv1[1:2]-mv2[1:2])
end

function rangeCompAllPoses(fgl1::AbstractDFG, fgl2::AbstractDFG)
  ranges = Float64[]
  xx,ll = ls(fgl1)
  for x in xx
    push!(ranges, rangeErrMaxPoint2(fgl1,x,fgl2,x))
  end
  return ranges
end

function rangeCompAllPoses(
    valsbaseline::Dict{Int,Array{Float64,1}},
    fglbaseline::AbstractDFG,
    fgltest::AbstractDFG)

  ranges = Float64[]
  xx,ll = ls(fgltest)
  for x in xx
    mv1 = valsbaseline[fglbaseline.IDs[x]]
    mv2 = getKDEMax(getVertKDE(fgltest,x))
    push!(ranges, norm(mv1[1:2]-mv2[1:2]))
  end
  return ranges
end





function selectFactorType(T1::Type{<:FunctorInferenceType}, T2::Type{<:FunctorInferenceType})
  # initial hacky version
  if T1 == Pose2 && T2 == Pose2
    return Pose2Pose2
  elseif T1 == Pose2 && T2 == Point2
    return Pose2Point2BearingRange
  elseif T1 == Point2 && T2 == Point2
    return Point2Point2
  else
    error("dont know which Factor type to select between $T1 and $T2")
  end
end

buildFactorDefault(::Type{Pose2Pose2}) = Pose2Pose2(MvNormal(zeros(3), diagm([0.01;0.01;0.01])))
buildFactorDefault(::Type{Point2Point2}) = Point2Point2(MvNormal(zeros(2), diagm([0.01;0.01])))
buildFactorDefault(::Type{Pose2Point2BearingRange}) = Pose2Point2BearingRange(Normal(0,0.1),Normal(1,0.1))
