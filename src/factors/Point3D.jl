

"""
$(TYPEDEF)

Direction observation information of a `Point3` variable.
"""
Base.@kwdef struct PriorPoint3{T} <: IncrementalInference.AbstractPrior where {T <: IIF.SamplableBelief}
  Z::T = MvNormal(zeros(3),diagm([1;1;1.0]))
end

getManifold(::InstanceType{PriorPoint3}) = getManifold(Point3)

# PriorPoint3 aka PriorPose3XYZ
function (cf::CalcFactor{<:PriorPoint3})(meas, X1::ArrayPartition)	
  return meas - X1.x[1]	
end

function (cf::CalcFactor{<:PriorPoint3})(meas, X1::AbstractVector)
  # return meas[1:3] .- X1[1:3] 	
  return meas - X1
end


## Serialized packing type and converters

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPriorPoint3  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end


function convert(::Type{PriorPoint3}, d::PackedPriorPoint3)
  return PriorPoint3(convert(SamplableBelief, d.Z))
end
function convert(::Type{PackedPriorPoint3}, d::PriorPoint3)
  return PackedPriorPoint3(convert(PackedSamplableBelief, d.Z))
end
