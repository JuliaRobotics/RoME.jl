
"""
    $TYPEDEF

Range and theta definition on `(:Euclid, :Circular)` manifold.
"""
struct Polar <: IIF.InferenceVariable
  dims::Int
  manifolds::Tuple{Symbol,Symbol}
  Polar() = new(2,(:Euclid,:Circular))
end

"""
    $TYPEDEF

Prior belief on any Polar related variable.
"""
mutable struct PriorPolar{T1<:IIF.SamplableBelief, T2<:IIF.SamplableBelief} <: IIF.AbstractPrior
  Zrange::T1
  Zangle::T2
  # empty constructor
  PriorPolar{T1,T2}() where {T1,T2} = new{T1,T2}()
  # regular constructor
  PriorPolar{T1,T2}(zr::T1, za::T2) where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = new{T1,T2}(zr,za)
end
# convenience and default object helper
PriorPolar(zr::T1=Normal(1,1), za::T2=Normal(0,0.1)) where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = PriorPolar{T1,T2}(zr,za)

function getSample(pp2i::PriorPolar, N::Int=1)
  sps = zeros(2,N)
  sps[1,:] = rand(pp2i.Zrange,N);
  sps[2,:] = rand(pp2i.Zangle,N);
  return (sps, )
end

"""
    $TYPEDEF

Linear offset factor of `IIF.SamplableBelief` between two `Polar` variables.
"""
mutable struct PolarPolar{T1<:IIF.SamplableBelief, T2<:IIF.SamplableBelief} <: IIF.AbstractRelativeRoots
  Zrange::T1
  Zangle::T2
  # empty constructor
  PolarPolar{T1,T2}() where {T1,T2} = new{T1,T2}()
  # regular constructor
  PolarPolar{T1,T2}(zr::T1, za::T2) where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = new{T1,T2}(zr,za)
end
# convenience and default object helper
PolarPolar(zr::T1=Normal(1,1), za::T2=Normal(0,0.1)) where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = PolarPolar{T1,T2}(zr,za)

function getSample(pp2i::PolarPolar, N::Int=1)
  sps = zeros(2,N)
  sps[1,:] = rand(pp2i.Zrange,N);
  sps[2,:] = rand(pp2i.Zangle,N);
  return (sps, )
end

function (pp::PolarPolar)(res::Array{Float64},
                          userdata::FactorMetadata,
                          idx::Int,
                          meas::Tuple{Array{Float64,2}},
                          p1::Array{Float64,2},
                          p2::Array{Float64,2} )
  #
  @inbounds res[1:2] = meas[1][1:2,idx] - (p2[1:2,idx] - p1[1:2,idx])
  nothing
end
