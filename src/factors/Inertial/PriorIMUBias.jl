
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorIMUBias( MvNormal(zeros(6), Matrix(Diagonal(ones(6).^2))) )
```
"""
Base.@kwdef struct PriorIMUBias{T <: SamplableBelief} <: IncrementalInference.AbstractPrior
  Z::T = MvNormal(zeros(6), diagm(0.5.*ones(6)))
end

DistributedFactorGraphs.getManifold(::InstanceType{PriorIMUBias}) = Manifolds.ProductGroup(
  ProductManifold(
    TranslationGroup(3), 
    TranslationGroup(3)
  )
)


# TODO the log here looks wrong (for gradients), consider:
# X = log(p⁻¹ ∘ m) 
# X = log(M, ϵ, Manifolds.compose(M, inv(M, p), m))
function (
  cf::CalcFactor{<:PriorIMUBias})(
  m::ArrayPartition, 
  p
)
  M = getManifold(PriorIMUBias)
  # TODO, Lie Group for now, expand to Riemannian
  ε = getPointIdentity(M)
  Xc = vee(M, ε, log(M, p, m))
  return Xc
end

#TODO Serialization of reference point p 
## Serialization support

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPriorIMUBias  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{PackedPriorIMUBias}, d::PriorIMUBias)
  return PackedPriorIMUBias(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorIMUBias}, d::PackedPriorIMUBias)
  return PriorIMUBias(convert(SamplableBelief, d.Z))
end

