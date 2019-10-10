# Graph analysis tools


function rangeErrMaxPoint2(fgl1::FactorGraph, id1, fgl2::FactorGraph ,id2)
  mv1 = getKDEMax(getVertKDE(fgl1,id1))
  mv2 = getKDEMax(getVertKDE(fgl2,id2))
  return norm(mv1[1:2]-mv2[1:2])
end

function rangeCompAllPoses(fgl1::FactorGraph, fgl2::FactorGraph)
  ranges = Float64[]
  xx,ll = ls(fgl1)
  for x in xx
    push!(ranges, rangeErrMaxPoint2(fgl1,x,fgl2,x))
  end
  return ranges
end

function rangeCompAllPoses(
    valsbaseline::Dict{Int,Array{Float64,1}},
    fglbaseline::FactorGraph,
    fgltest::FactorGraph)

  ranges = Float64[]
  xx,ll = ls(fgltest)
  for x in xx
    mv1 = valsbaseline[fglbaseline.IDs[x]]
    mv2 = getKDEMax(getVertKDE(fgltest,x))
    push!(ranges, norm(mv1[1:2]-mv2[1:2]))
  end
  return ranges
end
