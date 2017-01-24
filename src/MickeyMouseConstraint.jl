# Only 2D for now, 3 features from two poses
# will add multi-modality soon

# using Distributions, TransformUtils

function getUvecScaleFeature2D(wTb, bTc, wTo)
  cTo = (wTb * bTc) \ wTo
  r = cTo[1:2,3]
  α = norm(r)
  return r, r/α, α
end

function getUvecScaleBaseline2D(wTbi, wTbj, bTc)
  ciTcj = (wTbi * bTc) \ (wTbj * bTc)
  b = ciTcj[1:2,3]
  β = norm(b)
  return b, b/β, β
end

# function getEstFeatDirectii(wTbi, wTobjs::Vector, bTc)
#   xirvecs = Vector{Tuple}()
#   for wTo in wTobjs
#     push!(xirvecs, getUvecScaleFeature2D(wTbi, bTc, wTo))
#   end
#   return xirvecs
# end

function residruvec()
  nothing
end

# two sets of three feature sightings and a body to camera lever arm transform
type MickeyMouse2D <: IncrementalInference.FunctorPairwise
  # treat these as angles
  xir1::Normal
  xir2::Normal
  xir3::Normal
  xjr1::Normal
  xjr2::Normal
  xjr3::Normal
  bTc
  zDim::Tuple{Int,Int,Int,Int,Int}
  # Takes five
  MickeyMouse2D(a,b,c,d,e,f,g) = new(a,b,c,d,e,f,g, (3,3,2,2,2) )
end


function getSample(mm2d::MickeyMouse2D, N::Int=1)
  meas = Array{Float64,3}(3,N,2)
  meas[1,:,1] = rand(mm2d.xir1,N)
  meas[2,:,1] = rand(mm2d.xir2,N)
  meas[3,:,1] = rand(mm2d.xir3,N)
  meas[1,:,2] = rand(mm2d.xjr1,N)
  meas[2,:,2] = rand(mm2d.xjr2,N)
  meas[3,:,2] = rand(mm2d.xjr3,N)
  return (meas,)
end

function (mm2d::MickeyMouse2D)(res::Array{Float64}, idx::Int, meas::Tuple,
          wAbi::Array{Float64,2},
          wAbj::Array{Float64,2},
          wAo1::Array{Float64,2},
          wAo2::Array{Float64,2},
          wAo3::Array{Float64,2}  )
  #
  wTbi, wTbj = SE2(wAbi[:,idx]), SE2(wAbj[:,idx])
  wTo = ( SE2([wAo1[:,idx];0.0]),SE2([wAo2[:,idx];0.0]),SE2([wAo3[:,idx];0.0]) )

  b, bhat, β = getUvecScaleBaseline2D(wTbi, wTbj, mm2d.bTc)

  rhat, resid = zeros(2), zeros(2)
  # res[:] = -1e-5
  res[3] = -3*β

  for i in 1:3
    # sightings from first pose
    r1, rhathat1, α1 = getUvecScaleFeature2D(wTbi, mm2d.bTc, wTo[i])
    rhat[1] = cos(meas[1][i,idx,1])
    rhat[2] = sin(meas[1][i,idx,1])
    resid[1:2] = rhat-rhathat1
    res[1] += dot(resid, resid)

    # sightingings from second pose
    r2, rhathat2, α2 = getUvecScaleFeature2D(wTbj, mm2d.bTc, wTo[i])
    rhat[1] = cos(meas[1][i,idx,2])
    rhat[2] = sin(meas[1][i,idx,2])
    resid[1:2] = rhat-rhathat2
    res[2] += dot(resid, resid)

    # baseline constraint
    res[3] += dot(r1, bhat)
    res[3] += dot(r2, -bhat)
  end
  nothing
end















#
