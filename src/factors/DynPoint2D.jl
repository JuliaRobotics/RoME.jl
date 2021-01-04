# 2D SLAM with velocity states


"""
$(TYPEDEF)
"""
mutable struct DynPoint2VelocityPrior{T <: SamplableBelief} <: AbstractPrior
  z::T
  DynPoint2VelocityPrior{T}() where {T <: SamplableBelief} = new{T}()
  DynPoint2VelocityPrior{T}(z1::T) where {T <: SamplableBelief} = new{T}(z1)
end
DynPoint2VelocityPrior(z1::T) where {T <: SamplableBelief} = DynPoint2VelocityPrior{T}(z1)
getSample(dp2v::DynPoint2VelocityPrior, N::Int=1) = (rand(dp2v.z,N), )


"""
$(TYPEDEF)
"""
mutable struct DynPoint2DynPoint2{T <: SamplableBelief} <: AbstractRelativeRoots
  z::T
  DynPoint2DynPoint2{T}() where {T <: SamplableBelief} = new{T}()
  DynPoint2DynPoint2(z1::T) where {T <: SamplableBelief} = new{T}(z1)
end
getSample(dp2dp2::DynPoint2DynPoint2, N::Int=1) = (rand(dp2dp2.z,N), )
function (dp2dp2::DynPoint2DynPoint2)(
            res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple,
            Xi::Array{Float64,2},
            Xj::Array{Float64,2}  )
  #
  z = meas[1][:,idx]
  xi, xj = Xi[:,idx], Xj[:,idx]
  dt = Dates.value(userdata.fullvariables[2].nstime - userdata.fullvariables[1].nstime)*1e-9   # roughly the intended use of userdata
  res[1:2] = z[1:2] - (xj[1:2] - (xi[1:2]+dt*xi[3:4]))
  res[3:4] = z[3:4] - (xj[3:4] - xi[3:4])
  nothing
end




"""
$(TYPEDEF)
"""
mutable struct Point2Point2Velocity{T <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
  z::T
  Point2Point2Velocity{T}() where {T <: IIF.SamplableBelief} = new{T}()
  Point2Point2Velocity(z1::T) where {T <: IIF.SamplableBelief} = new{T}(z1)
end
getSample(p2p2v::Point2Point2Velocity, N::Int=1) = (rand(p2p2v.z,N), )
function (p2p2v::Point2Point2Velocity)(
                res::Array{Float64},
                userdata,
                idx::Int,
                meas::Tuple,
                Xi::Array{Float64,2},
                Xj::Array{Float64,2}  )
  #
  z = meas[1][:,idx]
  xi, xj = Xi[:,idx], Xj[:,idx]
  dt = (userdata.fullvariables[2].nstime - userdata.fullvariables[1].nstime)*1e-9     # roughly the intended use of userdata
  dp = (xj[1:2]-xi[1:2])
  dv = (xj[3:4]-xi[3:4])
  res[1] = 0.0
  res[1] += sum((z[1:2] - dp).^2)
  res[1] += sum((dp/dt - 0.5*(xj[3:4]+xi[3:4])).^2)  # (dp/dt - 0.5*(xj[3:4]+xi[3:4])) # midpoint integration
  res[1]
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
  return PackedDynPoint2VelocityPrior(string(d.z))
end
function convert(::Type{DynPoint2VelocityPrior}, d::PackedDynPoint2VelocityPrior)
  distr = extractdistribution(d.str)
  return DynPoint2VelocityPrior(distr)
end



#
