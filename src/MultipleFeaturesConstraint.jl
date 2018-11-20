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

"""
$(TYPEDEF)

two sets of three feature sightings and a body to camera lever arm transform
"""
mutable struct MultipleFeatures2D <: IncrementalInference.FunctorPairwiseMinimize
  # treat these as angles positive X->Y, X forward, Y left and Z up
  xir1::Normal
  xir2::Normal
  xir3::Normal
  xjr1::Normal
  xjr2::Normal
  xjr3a::Normal
  xjr3b::Normal
  abC::Categorical
  bTc::Array{Float64,2}
  zDim::Tuple{Int,Int,Int,Int,Int}
  # Takes five
  MultipleFeatures2D(a,b,c,d,e,f,g,h,i ) = new(a,b,c,d,e,f,g,h,i, (3,3,2,2,2) )
end

function getSample(mm2d::MultipleFeatures2D, N::Int=1)

  meas = Array{Float64,3}(undef, 3,N,2)
  meas[1,:,1] = rand(mm2d.xir1,N)
  meas[2,:,1] = rand(mm2d.xir2,N)
  meas[3,:,1] = rand(mm2d.xir3,N)
  meas[1,:,2] = rand(mm2d.xjr1,N)
  meas[2,:,2] = rand(mm2d.xjr2,N)
  meas[3,:,2] = rand(mm2d.xjr3a,N)

  select = rand(mm2d.abC,N)
  beb = select .== 2
  numb = sum(beb)
  meas[3,beb,2] = rand(mm2d.xjr3b,numb)

  return (meas, beb, numb)
end




# redo with angles on sightings for minimization
# NOTE -- an inefficient implementation
function (mm2d::MultipleFeatures2D)(res::Array{Float64},
          userdata,
          idx::Int, meas::Tuple,
          wAbi::Array{Float64,2},
          wAbj::Array{Float64,2},
          wAo1::Array{Float64,2},
          wAo2::Array{Float64,2},
          wAo3::Array{Float64,2},
          wAo3b::Union{Array{Float64,2},Bool}=false  )
  #
  len3b = wAo3b==false ? 0 : size(wAo3b,2)
  dobimodal = len3b > 0
  meas[3] > 0 && !dobimodal ? error("mm2d::MultipleFeatures2D -- bi-modal mismatch") : nothing

  wTbi = SE2(wAbi[:,idx])
  wTbj = SE2(wAbj[:,idx])
  wTo1 = SE2([wAo1[:,idx];0.0])
  wTo2 = SE2([wAo2[:,idx];0.0])
  wTo3 = SE2([wAo3[:,idx];0.0])
  wTo3b = dobimodal ? SE2([wAo3b[:,idx];0.0]) : nothing

  ciTo1 = (wTbi * mm2d.bTc) \ wTo1
  ciTo2 = (wTbi * mm2d.bTc) \ wTo2
  ciTo3 = (wTbi * mm2d.bTc) \ wTo3

  cjTo1 = (wTbj * mm2d.bTc) \ wTo1
  cjTo2 = (wTbj * mm2d.bTc) \ wTo2
  cjTo3 = (wTbj * mm2d.bTc) \ wTo3
  cjTo3b = dobimodal ? (wTbj * mm2d.bTc) \ wTo3 : nothing

  the1 = atan(se2vee(ciTo1)[2], se2vee(ciTo1)[1])
  the2 = atan(se2vee(ciTo2)[2], se2vee(ciTo2)[1])
  the3 = atan(se2vee(ciTo3)[2], se2vee(ciTo3)[1])

  res[1] = 0.0
  res[1] += ( meas[1][1,idx,1] - the1 )^2
  res[1] += ( meas[1][2,idx,1] - the2 )^2
  res[1] += ( meas[1][3,idx,1] - the3 )^2

  the1 = atan(se2vee(cjTo1)[2], se2vee(cjTo1)[1])
  the2 = atan(se2vee(cjTo2)[2], se2vee(cjTo2)[1])
  the3 = atan(se2vee(cjTo3)[2], se2vee(cjTo3)[1])
  the3b = dobimodal ? atan(se2vee(cjTo3b)[2], se2vee(cjTo3b)[1]) : nothing

  # res[2] = 0.0
  tempres = 0.0
  tempres += ( meas[1][1,idx,2] - the1 )^2
  tempres += ( meas[1][2,idx,2] - the2 )^2

  if dobimodal
    if !meas[2][idx]
      # 3a case
      tempres += ( meas[1][3,idx,2] - the3 )^2
    else
      # 3b case
      tempres += ( meas[1][3,idx,2] - the3b )^2
    end
  else
    tempres += ( meas[1][3,idx,2] - the3 )^2
  end
  res[1]+tempres
end
















#
