# 2D SLAM with velocity states

# x, y, dx/dt, dy/dt
mutable struct DynPoint2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  labels::Vector{String}
  DynPoint2(;ut::Int64=0, labels::Vector{<:AbstractString}=String[]) = new(ut, 4, labels)
end



mutable struct DynPoint2DynPoint2{T} <: IncrementalInference.FunctorPairwise where {T <: Distribution}
  z::T
  DynPoint2DynPoint2{T}() where {T <: Distribution} = new{T}()
  DynPoint2DynPoint2(z1::T) where {T <: Distribution} = new{T}(z1)
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
  dt = (userdata.variableuserdata[2].ut - userdata.variableuserdata[1].ut)*1e-6   # roughly the intended use of userdata
  res[1:2] = z[1:2] - (xj[1:2] - (xi[1:2]+dt*xi[3:4]))
  res[3:4] = z[3:4] - (xj[3:4] - xi[3:4])
  nothing
end
function (dp2dp2::DynPoint2DynPoint2)(res::Array{Float64},
            idx::Int,
            meas::Tuple,
            Xi::Array{Float64,2},
            Xj::Array{Float64,2}  )
  #
  error("function (dp2dp2::DynPoint2DynPoint2) requires newer version of IncrementalInference that passes user data to the residual function.")
end


mutable struct VelPoint2VelPoint2{T} <: IncrementalInference.FunctorPairwise where {T <: Distribution}
  z::T
  VelPoint2VelPoint2{T}() where {T <: Distribution} = new{T}()
  VelPoint2VelPoint2(z1::T) where {T <: Distribution} = new{T}(z1)
end
getSample(vp2vp2::VelPoint2VelPoint2, N::Int=1) = (rand(vp2vp2.z,N), )
function (vp2vp2::VelPoint2VelPoint2)(
                res::Array{Float64},
                userdata,
                idx::Int,
                meas::Tuple,
                Xi::Array{Float64,2},
                Xj::Array{Float64,2}  )
  #
  z = meas[1][:,idx]
  xi, xj = Xi[:,idx], Xj[:,idx]
  dt = (userdata.variableuserdata[2].ut - userdata.variableuserdata[1].ut)*1e-6   # roughly the intended use of userdata
  dp = (xj[1:2]-xi[1:2])
  res[1:2] = z[1:2] - dp
  res[3:4] = z[3:4] - (dp/dt - xi[3:4])
  # res[3:4] = z[3:4] - (dp/dt - 0.5*(xj[3:4]+xi[3:4])) # midpoint integration
  nothing # don't want to return anything from this function
end
function (vp2vp2::VelPoint2VelPoint2)(
                res::Array{Float64},
                idx::Int,
                meas::Tuple,
                Xi::Array{Float64,2},
                Xj::Array{Float64,2}  )
  #
    error("function (vp2vp2::VelPoint2VelPoint2) requires newer version of IncrementalInference that passes user data to the residual function.")
end


mutable struct DynPoint2VelocityPrior{T} <: IncrementalInference.FunctorSingleton where {T <: Distribution}
  z::T
  DynPoint2VelocityPrior{T}() where {T <: Distribution} = new{T}()
  DynPoint2VelocityPrior(z1::T) where {T <: Distribution} = new{T}(z1)
end
getSample(dp2v::DynPoint2VelocityPrior, N::Int=1) = (rand(dp2v.z,N), )
