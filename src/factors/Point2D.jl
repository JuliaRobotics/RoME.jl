
"""
$(TYPEDEF)

Direction observation information of a `Point2` variable.
"""
Base.@kwdef struct PriorPoint2{T <: IIF.SamplableBelief} <: IncrementalInference.AbstractPrior
  Z::T = MvNormal(zeros(2),LinearAlgebra.diagm([0.01;0.01]))
end

DFG.getManifold(::InstanceType{PriorPoint2}) = getManifold(Point2)


function (cf::CalcFactor{<:PriorPoint2})(meas, 	
                                        X1)	
  #	
  return meas[1:2] .- X1[1:2] 	
end

"""
$(TYPEDEF)
"""
Base.@kwdef struct Point2Point2{D <: IIF.SamplableBelief} <: AbstractManifoldMinimize #RelativeRoots
    Z::D = MvNormal(zeros(2),LinearAlgebra.diagm([0.1;0.1]))
end

getManifold(::InstanceType{Point2Point2}) = getManifold(Point2)


function (pp2r::CalcFactor{<:Point2Point2})(meas,
                                            xi,
                                            xj )
  #
  return meas[1:2] .- (xj[1:2] .- xi[1:2])
end




# ---------------------------------------------------------



"""
$(TYPEDEF)

Serialization type for `PriorPoint2`.
"""
Base.@kwdef struct PackedPriorPoint2  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end

function convert(::Type{PriorPoint2}, d::PackedPriorPoint2)
  return PriorPoint2(convert(SamplableBelief, d.Z))
end
function convert(::Type{PackedPriorPoint2}, d::PriorPoint2)
  return PackedPriorPoint2(convert(PackedSamplableBelief, d.Z))
end





"""
$(TYPEDEF)

Serialization type for `Point2Point2`.
"""
Base.@kwdef struct PackedPoint2Point2 <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{Point2Point2}, d::PackedPoint2Point2)
  return Point2Point2( convert(SamplableBelief, d.Z) )
end
function convert(::Type{PackedPoint2Point2}, d::Point2Point2)
  return PackedPoint2Point2( convert(PackedSamplableBelief, d.Z) )
end

