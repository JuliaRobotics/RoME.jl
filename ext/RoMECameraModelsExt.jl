module RoMECameraModelsExt

@info "RoME.jl is loading extension functionality using CameraModels.jl"

using CameraModels

import IncrementalInference: CalcFactor
import RoME: GenericProjection


function (CalcFactor{<:GenericProjection{S,T}})(
  c_X,
  w_P_c,
  w_P_o
) where {S,T}
  # predicted projection
  c_Xhat = project(S, T, w_P_c, w_P_o)
  # error vs measured projection
  return c_X - c_Xhat
end

end