

"""
$(TYPEDEF)
"""
Base.@kwdef struct Point3Point3{D <: IIF.SamplableBelief} <: AbstractRelativeMinimize
    Z::D = MvNormal(zeros(3),diagm([0.1;0.1;0.1]))
end

getManifold(::InstanceType{Point3Point3}) = getManifold(Point3)

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
Base.@kwdef struct PackedPoint3Point3 <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{Point3Point3}, d::PackedPoint3Point3)
  return Point3Point3( convert(SamplableBelief, d.Z) )
end
function convert(::Type{PackedPoint3Point3}, d::Point3Point3)
  return PackedPoint3Point3( convert(PackedSamplableBelief, d.Z) )
end