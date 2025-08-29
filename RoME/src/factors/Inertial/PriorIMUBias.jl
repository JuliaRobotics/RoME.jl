
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorIMUBias( MvNormal(zeros(6), Matrix(Diagonal(ones(6).^2))) )
```
"""
DFG.@defObservationType PriorIMUBias PriorObservation TranslationGroup(3) × TranslationGroup(3)

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
  Xc = m .- p
  return Xc
end

#TODO Serialization of reference point p 
## Serialization support


