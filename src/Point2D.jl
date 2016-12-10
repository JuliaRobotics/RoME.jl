# Point2D type


type PriorPoint2D <: IncrementalInference.Singleton
  mv::MvNormal
  W::Array{Float64,1}
  PriorPoint2D() = new()
  PriorPoint2D(mu, cov, W) = new(MvNormal(mu, cov), W)
end

type Point2DPoint2DRange <: IncrementalInference.Pairwise
    Zij::Vector{Float64} # bearing and range hypotheses as columns
    Cov::Float64
    W::Vector{Float64}
    Point2DPoint2DRange() = new()
    Point2DPoint2DRange(x...) = new(x[1],x[2],x[3])
end



function evalPotential(prior::PriorPoint2D, Xi::Array{Graphs.ExVertex,1}; N::Int64=100)#, from::Int64)
    return rand(prior.mv, N)
end


# Solve for Xid, given values from vertices [Xi] and measurement rho
function evalPotential(rho::Point2DPoint2DRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
  fromX, ret = nothing, nothing
  if Xi[1].index == Xid
    fromX = getVal( Xi[2] )
    ret = deepcopy(getVal( Xi[1] )) # carry pose yaw row over if required
  elseif Xi[2].index == Xid
    fromX = getVal( Xi[1] )
    ret = deepcopy(getVal( Xi[2] )) # carry pose yaw row over if required
  end
  r,c = size(fromX)
  theta = 2*pi*rand(c)
  noisy = rho.Cov*randn(c) + rho.Zij[1]

  for i in 1:c
    ret[1,i] = noisy[i]*cos(theta[i]) + fromX[1,i]
    ret[2,i] = noisy[i]*sin(theta[i]) + fromX[2,i]
  end
  return ret
end

























# ---------------------------------------------------------



type PackedPriorPoint2D  <: IncrementalInference.PackedInferenceType
    mu::Array{Float64,1}
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedPriorPoint2D() = new()
    PackedPriorPoint2D(x...) = new(x[1], x[2], x[3], x[4])
end

passTypeThrough(d::FunctionNodeData{Point2DPoint2DRange}) = d

function convert(::Type{PriorPoint2D}, d::PackedPriorPoint2D)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPoint2D(d.mu, Cov, d.W)
end
function convert(::Type{PackedPriorPoint2D}, d::PriorPoint2D)
  v2 = d.mv.Σ.mat[:];
  return PackedPriorPoint2D(d.mv.μ, v2, size(d.mv.Σ.mat,1), d.W)
end

function convert(::Type{PackedFunctionNodeData{PackedPriorPoint2D}}, d::FunctionNodeData{PriorPoint2D})
  return PackedFunctionNodeData{PackedPriorPoint2D}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          string(d.frommodule), convert(PackedPriorPoint2D, d.fnc))
end
function convert(::Type{FunctionNodeData{PriorPoint2D}}, d::PackedFunctionNodeData{PackedPriorPoint2D})
  return FunctionNodeData{PriorPoint2D}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
          Symbol(d.frommodule), convert(PriorPoint2D, d.fnc))
end
function FNDencode(d::FunctionNodeData{PriorPoint2D})
  return convert(PackedFunctionNodeData{PackedPriorPoint2D}, d)
end
function FNDdecode(d::PackedFunctionNodeData{PackedPriorPoint2D})
  return convert(FunctionNodeData{PriorPoint2D}, d)
end
