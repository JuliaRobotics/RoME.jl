
"""
$(TYPEDEF)

Direction observation information of a `Point2` variable.
"""
mutable struct PriorPoint2{T} <: IncrementalInference.AbstractPrior where {T <: IIF.SamplableBelief}
  Z::T
  # empty constructor
  PriorPoint2{T}() where T = new()
  # regular constructor
  PriorPoint2{T}(dist::T) where {T <: IIF.SamplableBelief} = new{T}(dist)
end
# convenience helper and default object
PriorPoint2(z::T=MvNormal(zeros(2),LinearAlgebra.diagm([0.01;0.01]))) where {T <: IIF.SamplableBelief} = PriorPoint2{T}(z)

function getSample(p2::PriorPoint2, N::Int=1)
  return (rand(p2.Z, N),)
end

#TODO wrapper
function (s::PriorPoint2{<:MvNormal})(X1::AbstractVector{T}; kwargs...) where T <: Real

  meas = mean(s.Z)
  iΣ = invcov(s.Z)
  res = meas[1:2] .- X1[1:2]
  return res' * iΣ * res

end

"""
$(TYPEDEF)
"""
mutable struct Point2Point2{D <: IIF.SamplableBelief} <: AbstractRelativeFactor
    Zij::D
    # empty constructor
    Point2Point2{T}() where T = new{T}()
    # regular constructor
    Point2Point2{T}(x::T) where {T <: IIF.SamplableBelief} = new{T}(x)
end
# convenience and default object helper
Point2Point2(x::T=MvNormal(zeros(2),LinearAlgebra.diagm([0.1;0.1]))) where {T <: IIF.SamplableBelief} = Point2Point2{T}(x)

function getSample(pp2::Point2Point2, N::Int=1)
  return (rand(pp2.Zij,N),  )
end
function (pp2r::Point2Point2{T})(
            res::Array{Float64},
            userdata::FactorMetadata,
            idx::Int,
            meas::Tuple,
            xi::Array{Float64,2},
            xj::Array{Float64,2} ) where T
  #
  res[1]  = meas[1][1,idx] - (xj[1,idx] - xi[1,idx])
  res[2]  = meas[1][2,idx] - (xj[2,idx] - xi[2,idx])
  nothing
end

#TODO wrapper
function (s::Point2Point2{<:MvNormal})(X1::AbstractVector{T}, X2::AbstractVector{T}; kwargs...) where T <: Real

  meas = mean(s.Zij)
  iΣ = invcov(s.Zij)
  res = meas[1:2] .- (X2[1:2] .- X1[1:2])
  return res' * iΣ * res
end





# ---------------------------------------------------------



"""
$(TYPEDEF)

Serialization type for `PriorPoint2`.
"""
mutable struct PackedPriorPoint2  <: IncrementalInference.PackedInferenceType
    str::String
    PackedPriorPoint2() = new()
    PackedPriorPoint2(x::String) = new(x)
end


function convert(::Type{PriorPoint2}, d::PackedPriorPoint2)
  # Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  distr = extractdistribution(d.str)
  return PriorPoint2{typeof(distr)}(distr)
end
function convert(::Type{PackedPriorPoint2}, d::PriorPoint2)
  # v2 = d.mv.Σ.mat[:];
  return PackedPriorPoint2(string(d.Z))
end





"""
$(TYPEDEF)

Serialization type for `Point2Point2`.
"""
mutable struct PackedPoint2Point2 <: IncrementalInference.PackedInferenceType
    str::String
    PackedPoint2Point2() = new()
    PackedPoint2Point2(s::AS) where {AS <: AbstractString} = new(s)
end
function convert(::Type{Point2Point2}, d::PackedPoint2Point2)
  return Point2Point2( extractdistribution(d.str) )
end
function convert(::Type{PackedPoint2Point2}, d::Point2Point2)
  return PackedPoint2Point2( string(d.Zij) )
end
