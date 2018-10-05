# 2D descriptor based on relative constraints in a marginalized map
# can beexpanded to 3D later.



mutable struct TreeDescriptor
  description::Array{Float64,2}
  meanidx::Float64
  position::Array{Float64,1}
end

treelms = Dict{Int,TreeDescriptor}()

# get the mean trajector index from sighing to this landmark
function getLandmTrajIndxs(fgl::FactorGraph, lmlbl::String)
  vlm = getVert(fgl,lmlbl)
  psoi = String[]

  for f in out_neighbors(fgl.g, vlm)
    for v2 in out_neighbors(fgl.g, f)
      if v2.attributes["label"][1] == "x" # use only poses
        # TODO -- add find first test here to ensure we do not double add pose
        push!(psoi, v2.attributes["label"])
      end
    end
  end
  return psoi
end

# take mean of poses connected to this landmark's index numbers
function getLandmMeanTrajIndx(fgl::FactorGraph, lmlbl::String)
  psoi = getLandmTrajIndxs(fgl, lmlbl)
  idxs = Float64[]
  for ps in psoi
    push!(idxs, parse(Float64, ps[2:end]))
  end
  return Base.mean(idxs)
end

function getAllLandmMeanTrajIndx(fgl::FactorGraph, LBS::Array{String,1})
  MIDX = Float64[]
  for i in 1:length(LBS)
    push!(MIDX, getLandmMeanTrajIndx(fgl, LBS[i]) )
  end
  return MIDX
end

# find idxs near mean idxs itol and range dtol
function getAllLandmNearByIndx(LB, X, Y, MIDX, nearidx;
                              itol::Float64=5.0, dtol::Float64=20.0)
  #permidx = Int[]
  # feat id and descriptor
  nnh = Dict{Int, TreeDescriptor}()

  x = X[nearidx]; y = Y[nearidx];
  dtol2 = dtol^2
  for i in 1:length(LB)
    xx = X[i]; yy = Y[i];
    d2 = (x-xx)^2 + (y-yy)^2
    if d2 < dtol2 && abs(MIDX[nearidx]-MIDX[i]) < itol && i != nearidx
      #push!(permidx, i)
      fid = parse(Int, LB[i][2:end])
      nnh[fid] = TreeDescriptor(Array{Float64,2}(), MIDX[i], [xx;yy])
    end
  end
  return nnh
end

function calcFeatDesc(fgl::FactorGraph, flbl::String, midx::Float64)
  lm = getVert(fgl,flbl)
end

function calcLandmDescriptions!(fgl::FactorGraph, lmd::Dict{Int,TreeDescriptor})
  X,Y,th,LB = get2DLandmMeans(fgl)
  getAllLandmNearByIndx(LB,X,Y)
  MIDX = getAllLandmMeanTrajIndx(fgl, LB)

  for i in 1:length(LB)
    featid = parse(Float64,LB[i][2:end])
    nnh = getAllLandmNearByIndx(LB, X, Y, MIDX, i) # this should be done at feature creation
    desc = calcFeatDesc(fgl, flbl, )
  end

end
