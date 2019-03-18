

"""
$(TYPEDEF)

Direction observation information of a `Point3` variable.
"""
mutable struct PriorPoint3{T} <: IncrementalInference.FunctorSingleton where {T <: IIF.SamplableBelief}
  Z::T
  # W::Array{Float64,1} # TODO, deprecate the weight parameter
  PriorPoint3{T}() where T = new()
  PriorPoint3{T}(dist::T) where {T <: IIF.SamplableBelief} = new{T}(dist)
end
PriorPoint3(z::T) where {T <: IIF.SamplableBelief} = PriorPoint3{T}(z)

function getSample(p3::PriorPoint3, N::Int=1)
  return (rand(p3.Z, N),)
end



## Serialized packing type and converters

"""
$(TYPEDEF)
"""
mutable struct PackedPriorPoint3  <: IncrementalInference.PackedInferenceType
    str::String
    PackedPriorPoint3() = new()
    PackedPriorPoint3(x::String) = new(x)
end


function convert(::Type{PriorPoint3}, d::PackedPriorPoint3)
  # Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  distr = extractdistribution(d.str)
  return PriorPoint3{typeof(distr)}(distr)
end
function convert(::Type{PackedPriorPoint3}, d::PriorPoint3)
  # v2 = d.mv.Σ.mat[:];
  return PackedPriorPoint3(string(d.Z))
end
