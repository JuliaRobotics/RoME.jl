
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
# mutable struct PriorPose2{T} <: IncrementalInference.AbstractPrior  where {T <: IncrementalInference.SamplableBelief}
#   Z::T
#   # empty constructor
#   PriorPose2{T}() where T = new{T}()
#   # regular constructor
#   PriorPose2{T}(x::T) where {T <: IncrementalInference.SamplableBelief}  = new{T}(x)
# end
# # convenience and default object helper
# PriorPose2(x::T) where {T <: IncrementalInference.SamplableBelief} = PriorPose2{T}(x)
struct PriorPose2{T <: SamplableBelief, P} <: IIF.AbstractPrior
  Z::T
  p::P 
end

PriorPose2(z::IIF.SamplableBelief) = PriorPose2(z, getPointIdentity(Pose2))
PriorPose2{T,P}(z::IIF.SamplableBelief) where {T,P} = PriorPose2{T,P}(z, getPointIdentity(Pose2))


DFG.getManifold(::PriorPose2) = SpecialEuclidean(2)


function getSample(cf::CalcFactor{<:PriorPose2}, N::Int=1)
  Z = cf.factor.Z
  p = cf.factor.p
  M = getManifold(cf.factor)
  
  Xc = rand(Z)
  
  # X = get_vector.(Ref(M), Ref(p), Xc, Ref(DefaultOrthogonalBasis()))
  X = hat(M, p, Xc)
  points = exp(M, p, X)

  return (points, )
end

function (cf::CalcFactor{<:PriorPose2})(m, p)	
  #	
  M = getManifold(cf.factor)
  return log(M, p, m)
  # iXihat = SE2(meas) \ SE2(wXi)	
  # return se2vee(iXihat)	
end

#TODO Serialization of reference point p 
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
  return PackedPriorPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  distr = convert(SamplableBelief, d.str)
  return PriorPose2(distr)
end




## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  compareDensity(a.Z, b.Z)
end
