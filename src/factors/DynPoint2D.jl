# 2D SLAM with velocity states


"""
$(TYPEDEF)
"""
mutable struct DynPoint2VelocityPrior{T <: SamplableBelief} <: AbstractPrior
  Z::T
end

getManifold(::DynPoint2VelocityPrior) = TranslationGroup(4)

"""
$(TYPEDEF)
"""
mutable struct DynPoint2DynPoint2{T <: SamplableBelief} <: AbstractRelativeRoots
  Z::T
end

getManifold(::DynPoint2DynPoint2) = TranslationGroup(4)


function (cfo::CalcFactor{<:DynPoint2DynPoint2})(z, xi, xj)
  #
  dt = Dates.value(cfo.fullvariables[2].nstime - cfo.fullvariables[1].nstime)*1e-9   # roughly the intended use of userdata
  res12 = z[1:2] - (xj[1:2] - (xi[1:2]+dt*xi[3:4]))
  res34 = z[3:4] - (xj[3:4] - xi[3:4])
  return [res12; res34]
end


"""
$(TYPEDEF)
"""
mutable struct Point2Point2Velocity{T <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
  Z::T
end

getManifold(::Point2Point2Velocity) = TranslationGroup(4)

function (cfo::CalcFactor{<:Point2Point2Velocity})( z,
                                                    xi,
                                                    xj  )
  #
  dt = (cfo.fullvariables[2].nstime - cfo.fullvariables[1].nstime)*1e-9     # roughly the intended use of userdata
  dp = (xj[1:2] .- xi[1:2])
  dv = (xj[3:4] .- xi[3:4])

  res12 = z[1:2] .- dp
  res34 =  dp/dt .- 0.5*(xj[3:4] .+ xi[3:4])  # (dp/dt - 0.5*(xj[3:4]+xi[3:4])) # midpoint integration

  return [res12; res34]
end



## Packing Types================================================================


"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedDynPoint2VelocityPrior <: AbstractPackedFactor
  str::PackedSamplableBelief
end

function convert(::Type{PackedDynPoint2VelocityPrior}, d::DynPoint2VelocityPrior)
  return PackedDynPoint2VelocityPrior(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{DynPoint2VelocityPrior}, d::PackedDynPoint2VelocityPrior)
  distr = convert(SamplableBelief, d.str)
  return DynPoint2VelocityPrior(distr)
end



#
