
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
mutable struct PriorPose2{T} <: IncrementalInference.AbstractPrior  where {T <: IncrementalInference.SamplableBelief}
  Z::T
  # empty constructor
  PriorPose2{T}() where T = new{T}()
  # regular constructor
  PriorPose2{T}(x::T) where {T <: IncrementalInference.SamplableBelief}  = new{T}(x)
end
# convencience and default object helper
PriorPose2(x::T) where {T <: IncrementalInference.SamplableBelief} = PriorPose2{T}(x)

function getSample(p2::PriorPose2, N::Int=1)
  return (rand(p2.Z,N), )
end

#TODO wrapper
function (s::PriorPose2{<:MvNormal})(wXi::AbstractVector{T}; kwargs...) where T <: Real


  meas = mean(s.Z)
  iΣ = invcov(s.Z)

  # res = meas .- X1
  iXihat = SE2(meas[1:3]) \ SE2(wXi[1:3])
  res = se2vee(iXihat)

  return res' * iΣ * res
end




## Serialization support

"""
$(TYPEDEF)
"""
mutable struct PackedPriorPose2  <: IncrementalInference.PackedInferenceType
    str::String
    PackedPriorPose2() = new()
    PackedPriorPose2(x::String) = new(x)
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  return PackedPriorPose2(string(d.Z))
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  distr = extractdistribution(d.str)
  return PriorPose2{typeof(distr)}(distr)
end




## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  TP = true
  TP = TP && norm(a.Z.μ-b.Z.μ) < tol
  TP = TP && norm(a.Z.Σ.mat-b.Z.Σ.mat) < tol
  return TP
end
