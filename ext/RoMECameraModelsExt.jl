module RoMECameraModelsExt

# @info "RoME.jl is loading extension functionality using CameraModels.jl"

using CameraModels
using StaticArrays
using Manifolds
using DocStringExtensions
using Optim

import IncrementalInference: AbstractDFG, getFactorType, getVariable, getSolverData, CalcFactor, ls
import RoME: GenericProjection, solveMultiviewLandmark!


function projectPointFrom(cam, c_H_w, w_Ph)
  c_Ph = c_H_w*w_Ph |> SVector{4}
  CameraModels.projectHomogeneous(cam,c_Ph)
end

function (cf::CalcFactor{<:GenericProjection{S,T}})(
  c_X,
  w_P_c,
  w_P_o
) where {S,T}
  
  κ = 1000
  M = SpecialEuclidean(3)

  w_Ph = SVector{4,Float64}(w_P_o..., 1.0)
  c_H_w = inv(affine_matrix(M, w_P_c))

  # predicted projection
  # c_Xhat = project(S, T, w_P_c, w_P_o)
  pred = projectPointFrom(cf.factor.cam, c_H_w, w_Ph)
  # experimental cost function to try force bad reprojects in front of the camera during optimization
  # κ*(abs(pred.depth) - pred.depth)^2 + (c_X[1]-pred[1])^2 + (c_X[2]-pred[2])^2

  res = SVector{3,Float64}(
    c_X[1]-pred[1],
    c_X[2]-pred[2],
    κ*(abs(pred.depth) - pred.depth)
  )

  # error vs measured projection
  # return c_X - c_Xhat
  return res
end







"""
    $SIGNATURES

native Optim solution to complement GenericProjection factor development.

Notes
- See RoME/test/testGenericProjection.jl for usage example.
- Original implementation from JuliaRobotics/CameraModels.jl/test/multiview_manifolds.jl
"""
function solveMultiviewLandmark!(
  dfg::AbstractDFG, 
  lmlb::Symbol;
  cam = CameraModels.CameraCalibration(),
  retry::Int = 100
)
  # assume to find 
  M = SpecialEuclidean(3)
  flbs = ls(dfg, lmlb)
  # fcs = getFactor.(dfg, flbs)
  # fcd = getFactorType.(fcs)
  
  function projectPointFrom(cam, c_H_w, w_Ph)
    c_Ph = c_H_w*w_Ph |> SVector{4}
    CameraModels.projectHomogeneous(cam,c_Ph)
  end
  
  function cameraResidual(cam, meas, M, w_T_c, w_Ph, κ=1000)
    pred = projectPointFrom(cam, inv(affine_matrix(M,w_T_c)), w_Ph)
    # experimental cost function to try force bad reprojects in front of the camera during optimization
    κ*(abs(pred.depth) - pred.depth)^2 + (meas[1]-pred[1])^2 + (meas[2]-pred[2])^2
  end
  
  function cost(cam, c_Xi, w_P_ci, w_Ph)
    res = 0.0
    for (c_X, w_P_c) in zip(c_Xi, w_P_ci)
      res += cameraResidual(cam, c_X, M, w_P_c, w_Ph)
    end
    return res
  end

  _w_P_ci = Vector{ArrayPartition{Float64, Tuple{SVector{3, Float64}, SMatrix{3, 3, Float64, 9}}}}()
  vlbs = Symbol[]
  _c_Xi = Vector{CameraModels.PixelIndex{true, Float64}}()
  for fl in flbs
    vl = setdiff(ls(dfg, fl), [lmlb;])[1]
    push!(
      _w_P_ci,
      getSolverData(getVariable(dfg, vl), :parametric).val[1]
    )
    union!(
      vlbs,
      [vl;],
    )
    
    obs = CameraModels.PixelIndex(getFactorType(dfg, fl).Z.μ...)
    push!(_c_Xi, obs)
  end
  cost_(w_Lh) = cost(cam, _c_Xi, _w_P_ci, SVector(w_Lh...))

  w_P3 = zeros(3)
  for i in 1:retry
    w_Ph0 = [(retry*randn(3) .+ getSolverData(getVariable(dfg, lmlb), :parametric).val[1])...; 1.0]
    # w_Ph0 = [10;0;0;1.0]
    w_Res = Optim.optimize(
      cost_, 
      w_Ph0, # [1;0;0;1.], 
      LBFGS(),
      # ParticleSwarm(; lower = [0.,-1.,-1.,0.],
      #               upper = [99999.;1;1;9999],
      #               n_particles = 10),
      # Optim.Options(g_tol = 1e-12,
      #               iterations = 100,
      #               store_trace = false,
      #               show_trace = true);
      # autodiff=:forward
    )
    w_P3 = w_Res.minimizer |> CameraModels.toNonhomogeneous
    # check depth okay
    w_Ph = [w_P3...; 1.0]
    depthok = true
    for w_T_c in _w_P_ci 
      pred = projectPointFrom(cam, inv(affine_matrix(M, w_T_c)), w_Ph)
      if pred.depth < 0
        depthok = false
      end
    end
    if w_Res.ls_success && depthok
      break
    end
    if retry <= i 
      throw(DomainError("Unable to converge projection solution"))
    end
  end

  return w_P3
end



end