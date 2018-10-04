
mutable struct CameraIntrinsic
  K::Array{Float64,2}
  CameraIntrinsic(::Nothing) = new()
  CameraIntrinsic(;x0=320.0,y0=240.0,fx=510.0,fy=510.0,s=0.0) = new([[fx;s;x0]';[0.0;fy;y0]';[0.0;0;1]'])
end

# Camera extrinsic must be world in camera frame (cRw)
mutable struct CameraExtrinsic
  R::SO3
  t::Vector{Float64}
  CameraExtrinsic(::Nothing) = new()
  CameraExtrinsic(;R=SO3(0),t=zeros(3)) = new(R, t)
end
mutable struct CameraModelFull
  ci::CameraIntrinsic
  ce::CameraExtrinsic
  # cd::CameraDistortion
  CameraModelFull(::Nothing) = new()
  CameraModelFull(;ci=CameraIntrinsic(), ce=CameraExtrinsic()) = new(ci,ce)
end
function project!(ret::Vector{Float64}, ci::CameraIntrinsic, ce::CameraExtrinsic, pt::Vector{Float64})
  res = ci.K*(ce.R.R*pt + ce.t)
  ret[1:2] = res[1:2]./res[3]
  nothing
end
project!(ret::Vector{Float64}, cm::CameraModelFull, pt::Vector{Float64}) = project!(ret, cm.ci, cm.ce, pt)
function project(cm::CameraModelFull, pt::Vector{Float64})
  res = Vector{Float64}(2)
  project!(res, cm, pt)
  return res
end

# pinhole camera model
# (x, y)/f = (X, Y)/Z
function cameraResidual!(
      res::Vector{Float64},
      z::Vector{Float64},
      ci::CameraIntrinsic,
      ce::CameraExtrinsic,
      pt::Vector{Float64}  )
  # in place memory operations
  project!(res, ci, ce, pt)
  res[1:2] .*= -1.0
  res[1:2] += z[1:2]
  nothing
end
