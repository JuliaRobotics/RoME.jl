
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
Base.@kwdef struct PriorPose2{T <: SamplableBelief, P} <: IIF.AbstractPrior
  Z::T
end

DFG.getManifold(::InstanceType{PriorPose2}) = getManifold(Pose2) # SpecialEuclidean(2)


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
Base.@kwdef struct PackedPriorPose2  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  return PackedPriorPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  return PriorPose2(convert(SamplableBelief, d.Z))
end




## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  compareDensity(a.Z, b.Z)
end
