










"""
$(TYPEDEF)
"""
mutable struct Point3Point3{D <: IIF.SamplableBelief} <: AbstractRelativeRoots
    Zij::D
    # empty constructor
    Point3Point3{T}() where T = new{T}()
    # regular constructor
    Point3Point3{T}(x::T) where {T <: IIF.SamplableBelief} = new{T}(x)
end
# convenience and default object helper
Point3Point3(x::T=MvNormal(zeros(3),LinearAlgebra.diagm([0.1;0.1;0.1]))) where {T <: IIF.SamplableBelief} = Point3Point3{T}(x)

function getSample(cfo::CalcFactor{<:Point3Point3}, N::Int=1)
  return ([rand(cfo.factor.Zij) for _=1:N],  )
end
function (cf::CalcFactor{<:Point3Point3})(meas,
                                            xi,
                                            xj )
  #
  return meas[1:3] .- (xj[1:3] .- xi[1:3])
end










"""
$(TYPEDEF)

Serialization type for `Point3Point3`.
"""
mutable struct PackedPoint3Point3 <: IncrementalInference.PackedInferenceType
    str::String
    PackedPoint3Point3() = new()
    PackedPoint3Point3(s::AS) where {AS <: AbstractString} = new(s)
end
function convert(::Type{Point3Point3}, d::PackedPoint3Point3)
  return Point3Point3( convert(SamplableBelief, d.str) )
end
function convert(::Type{PackedPoint3Point3}, d::Point3Point3)
  return PackedPoint3Point3( convert(PackedSamplableBelief, d.Zij) )
end