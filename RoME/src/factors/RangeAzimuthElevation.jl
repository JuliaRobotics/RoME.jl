
# FIXME update to Manifolds.Sphere instead

mutable struct RangeAzimuthElevation
  range::Float64
  azimuth::Float64
  elevation::Union{Nothing,Float64}
end


function convert(::Type{RangeAzimuthElevation}, val::Tuple{Symbol, Vector{Float64}})
  if val[1] == :rangeazimuth
    return RangeAzimuthElevation(val[2][1],val[2][2],nothing)
  elseif val[1] == :rangeazimuthelevation
    return RangeAzimuthElevation(val[2][1],val[2][2],val[2][3])
  else
    error("Unknown conversion from $(val[1]) to RangeAzimuthElevation.")
  end
end


function \(s::SE3, wTr::CTs.Translation)
  bTr = s.R.R'*(wTr.v-s.t)
  Dtr = bTr
  range = norm(Dtr)
  azi = atan(Dtr[2], Dtr[1])
  elev = atan(Dtr[3], Dtr[1])
  RangeAzimuthElevation(range, azi, elev)
end
