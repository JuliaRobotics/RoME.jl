# 2D SLAM with velocity states

# x, y, dx/dt, dy/dt
mutable struct DynPoint2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  labels::Vector{String}
  DynPoint2(;ut::Int64=0, labels::Vector{<:AbstractString}=String[]) = new(ut, 4, lbls)
end



mutable struct DynPoint2DynPoint2{T} <: IncrementalInference.FunctorPairwise where {T <: Distribution}
  z::T
  DynPoint2DynPoint2{T}() where {T <: Distribution} = new{T}()
  DynPoint2DynPoint2(z1::T) where {T <: Distribution} = new{T}(z1)
end
getSample(dp2dp2::DynPoint2DynPoint2, N::Int=1) = (rand(s.z,N), )
function (dp2dp2::DynPoint2DynPoint2)(
            res::Array{Float64},
            idx::Int,
            userdata ,
            meas::Tuple,
            Xi::Array{Float64,2},
            Xj::Array{Float64,2}  )
  #
  Z = meas[1][1,idx]
  xi, xj = Xj[:,idx],Xj[:,idx]
  # dt = xj[5] - xi[5]
  dt = NaN # dt = (userdata.tj - userdata.ti)*1e-6   # roughly the intended use of userdata
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
