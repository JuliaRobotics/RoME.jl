## ==================================================================================================
## Hack to be removed or updated
## FIXME ME ON FIRE use the ManifoldsBase.jl prescribed interface method instead:
##   https://juliamanifolds.github.io/Manifolds.jl/stable/examples/manifold.html#manifold-tutorial
## ==================================================================================================


import ApproxManifoldProducts: coords, uncoords, getPointsManifold, _makeVectorManifold

export SE2E2_Manifold


# this is a hack and not fully implemented as per: 
struct _SE2E2 <: ManifoldsBase.Manifold{ManifoldsBase.ℝ} end

const SE2E2_Manifold = _SE2E2()

ManifoldsBase.manifold_dimension(::_SE2E2) = 5

AMP.coords(::Type{<:typeof(SE2E2_Manifold)}, p::Manifolds.ProductRepr) = [p.parts[1][1], p.parts[1][2], atan(p.parts[2][2,1],p.parts[2][1,1]), p.parts[3][1], p.parts[3][2]]

function AMP.uncoords(::typeof(SE2E2_Manifold), p::AbstractVector{<:Real})
  α = p[3]
  Manifolds.ProductRepr(([p[1], p[2]]), [cos(α) -sin(α); sin(α) cos(α)], ([p[4], p[5]]))
end

function AMP.getPointsManifold(mkd::ManifoldKernelDensity{M}) where {M <: typeof(SE2E2_Manifold)}
  data_ = getPoints(mkd.belief)
  [uncoords(mkd.manifold, view(data_, :, i)) for i in 1:size(data_,2)]
end

function Statistics.mean(::typeof(SE2E2_Manifold), pts::AbstractVector)
  se2_ = (d->Manifolds.ProductRepr(d.parts[1], d.parts[2])).(pts)
  mse2 = mean(Manifolds.SpecialEuclidean(2), se2_)
  e2_ = (d->Manifolds.ProductRepr(d.parts[3])).(pts)
  me2 = mean(Euclidean(2), e2_)
  Manifolds.ProductRepr(mse2.parts[1], mse2.parts[2], me2.parts[1])
end

AMP._makeVectorManifold(::M, prr::ProductRepr) where {M <: typeof(SE2E2_Manifold)} = coords(M, prr)




## =============================================================================
## Needs a home
## =============================================================================

export BearingRange_Manifold


struct _CircleEuclid <: ManifoldsBase.Manifold{ManifoldsBase.ℝ} end

ManifoldsBase.manifold_dimension(::_CircleEuclid) = 2

const BearingRange_Manifold = _CircleEuclid()

AMP.coords(::Type{<:typeof(BearingRange_Manifold)}, p::ProductRepr) = [p.parts[1][1]; p.parts[2][1]]

function AMP.uncoords(::typeof(BearingRange_Manifold), p::AbstractVector{<:Real})
  ProductRepr(([p[1];]), ([p[2];]))
end

function AMP.getPointsManifold(mkd::ManifoldKernelDensity{M}) where {M <: typeof(BearingRange_Manifold)}
  data_ = getPoints(mkd.belief)
  [uncoords(mkd.manifold, view(data_, :, i)) for i in 1:size(data_,2)]
end

function Statistics.mean(::typeof(BearingRange_Manifold), pts::AbstractVector)
  TensorCast.@cast bearing[i] := pts[i][1]
  TensorCast.@cast range_[i] := pts[i][2]
  mc = mean(Circle(), bearing)
  mr = mean(range_)

  return [mc; mr]
end

AMP._makeVectorManifold(::M, prr::ProductRepr) where {M <: typeof(BearingRange_Manifold)} = coords(M, prr)



## ==================================================================================================
