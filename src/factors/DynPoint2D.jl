# 2D SLAM with velocity states


"""
$(TYPEDEF)
"""
mutable struct DynPoint2VelocityPrior{T <: SamplableBelief} <: AbstractPrior
  Z::T
  DynPoint2VelocityPrior{T}() where {T <: SamplableBelief} = new{T}()
  DynPoint2VelocityPrior{T}(z1::T) where {T <: SamplableBelief} = new{T}(z1)
end
DynPoint2VelocityPrior(z1::T) where {T <: SamplableBelief} = DynPoint2VelocityPrior{T}(z1)

getSample(cfo::CalcFactor{<:DynPoint2VelocityPrior}) = rand(cfo.factor.Z)


"""
$(TYPEDEF)
"""
mutable struct DynPoint2DynPoint2{T <: SamplableBelief} <: AbstractRelativeRoots
  Z::T
  DynPoint2DynPoint2{T}() where {T <: SamplableBelief} = new{T}()
  DynPoint2DynPoint2(z1::T) where {T <: SamplableBelief} = new{T}(z1)
end

getSample(cfo::CalcFactor{<:DynPoint2DynPoint2}) = rand(cfo.factor.Z)

function (cfo::CalcFactor{<:DynPoint2DynPoint2})(Z, xi, xj)
  #
  dt = Dates.value(cfo.metadata.fullvariables[2].nstime - cfo.metadata.fullvariables[1].nstime)*1e-9   # roughly the intended use of userdata
  res12 = Z[1:2] - (xj[1:2] - (xi[1:2]+dt*xi[3:4]))
  res34 = Z[3:4] - (xj[3:4] - xi[3:4])
  return [res12; res34]
end


"""
$(TYPEDEF)
"""
mutable struct Point2Point2Velocity{T <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
  Z::T
  Point2Point2Velocity{T}() where {T <: IIF.SamplableBelief} = new{T}()
  Point2Point2Velocity(z1::T) where {T <: IIF.SamplableBelief} = new{T}(z1)
end

getSample(cfo::CalcFactor{<:Point2Point2Velocity}) = rand(cfo.factor.Z)
function (cfo::CalcFactor{<:Point2Point2Velocity})( Z,
                                                    xi,
                                                    xj  )
  #
  dt = (cfo.metadata.fullvariables[2].nstime - cfo.metadata.fullvariables[1].nstime)*1e-9     # roughly the intended use of userdata
  dp = (xj[1:2] .- xi[1:2])
  dv = (xj[3:4] .- xi[3:4])

  res12 = Z[1:2] .- dp
  res34 =  dp/dt .- 0.5*(xj[3:4] .+ xi[3:4])  # (dp/dt - 0.5*(xj[3:4]+xi[3:4])) # midpoint integration

  return [res12; res34]
end



## Packing Types================================================================


"""
$(TYPEDEF)
"""
mutable struct PackedDynPoint2VelocityPrior <: IIF.PackedInferenceType
  str::String
  PackedDynPoint2VelocityPrior() = new()
  PackedDynPoint2VelocityPrior(z1::String) = new(z1)
end

function convert(::Type{PackedDynPoint2VelocityPrior}, d::DynPoint2VelocityPrior)
  return PackedDynPoint2VelocityPrior(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{DynPoint2VelocityPrior}, d::PackedDynPoint2VelocityPrior)
  distr = convert(SamplableBelief, d.str)
  return DynPoint2VelocityPrior(distr)
end



#
