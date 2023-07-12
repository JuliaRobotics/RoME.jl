
"""
    $TYPEDEF

Range and theta definition on `(:Euclid, :Circular)` manifold.
"""
struct Polar <: IIF.InferenceVariable
  dims::Int
end

getManifold(::InstanceType{Polar}) = BearingRange_Manifold

"""
    $TYPEDEF

Prior belief on any Polar related variable.
"""
Base.@kwdef struct PriorPolar{T1<:IIF.SamplableBelief, T2<:IIF.SamplableBelief} <: IIF.AbstractPrior
  Zrange::T1 = Normal(1,1)
  Zangle::T2 = Normal(0,0.1)
end

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
Base.@kwdef struct PolarPolar{T1<:IIF.SamplableBelief, T2<:IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
  Zrange::T1 = Normal(1,1)
  Zangle::T2 = Normal(0,0.1)
end

function getSample(pp2i::PolarPolar)
  sps = zeros(2)
  sps[1] = rand(pp2i.Zrange);
  sps[2] = rand(pp2i.Zangle);
  return sps
end

function (pp::PolarPolar)(meas,
                          p1,
                          p2 )
  #
  return meas .- (p2 .- p1)
end
