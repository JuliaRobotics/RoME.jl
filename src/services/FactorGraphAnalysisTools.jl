# Graph analysis tools and support calculation tools


function rangeErrMaxPoint2(fgl1::AbstractDFG, id1, fgl2::AbstractDFG ,id2)
  mv1 = getKDEMax(getBelief(fgl1,id1))
  mv2 = getKDEMax(getBelief(fgl2,id2))
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
    mv2 = getKDEMax(getBelief(fgltest,x))
    push!(ranges, norm(mv1[1:2]-mv2[1:2]))
  end
  return ranges
end


## Compare functions

function compareDensity(a::Normal, b::Normal; tol::Float64=1e-10)
  TP = true
  TP = TP && abs(a.μ - b.μ) < tol
  TP = TP && abs(a.σ - b.σ) < tol
  TP
end

function compareDensity(a::MvNormal, b::MvNormal; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && norm(a.μ - b.μ)<tol
  TP = TP && sum(norm.(a.Σ.mat - b.Σ.mat))<tol
  return TP
end


#
