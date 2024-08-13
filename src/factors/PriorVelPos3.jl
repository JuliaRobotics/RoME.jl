
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorVelPos3( MvNormal(zeros(6), Matrix(Diagonal(ones(6).^2))) )
```
"""
Base.@kwdef struct PriorVelPos3{T <: SamplableBelief} <: IIF.AbstractPrior
  Z::T = MvNormal(zeros(3), diagm([1;1;0.1]))
end

DFG.getManifold(::InstanceType{PriorVelPos3}) = getManifold(VelPos3)


# TODO the log here looks wrong (for gradients), consider:
# X = log(p⁻¹ ∘ m) 
# X = log(M, ϵ, Manifolds.compose(M, inv(M, p), m))
function (
  cf::CalcFactor{<:PriorVelPos3})(
  m::ArrayPartition, 
  p::ArrayPartition
)
  M = getManifold(PriorVelPos3)
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
Base.@kwdef struct PackedPriorVelPos3  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{PackedPriorVelPos3}, d::PriorVelPos3)
  return PackedPriorVelPos3(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorVelPos3}, d::PackedPriorVelPos3)
  return PriorVelPos3(convert(SamplableBelief, d.Z))
end

