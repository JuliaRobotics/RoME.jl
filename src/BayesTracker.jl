
# TODO -- make NBP a module
# include("../NBP/BallTree01.jl")
# include("../NBP/BallTreeDensity01.jl")
# include("../NBP/KDE01.jl")
# include("../NBP/DualTree01.jl")
# include("../NBP/MSGibbs01.jl")
# include("../NBP/KDEPlotting01.jl")

type Feature
  id::Int64
  age::Int64
  lastzage::Int64
  lastz::Array{Float64,1}
  bel::BallTreeDensity
end

locpid=1
function upplocal()
    global locpid
    N = nprocs()
    locpid = (N > 1 ? ((locpid < N && locpid != 0 && N > 2) ? (locpid+1)%(N+1) : 2) : (N == 1 ? 1 : 2))
    return locpid
end

# trackers = Dict{Int64,Feature}()


function propagate(b1Xl::BallTreeDensity, bDxb1::Array{Float64,1}, s::Array{Float64,1})
  pts1 = getPoints(b1Xl)
  pts = zeros(size(pts1))
  for i in 1:size(pts,2)
    ent = [s[1]*randn();s[2]*randn();s[3]*randn()]
    b1Tb = inv(SE2(bDxb1+ent))
    pt = vec(pts1[:,i])
    bTl = SE2([pt[1];pt[2];0.0])
    pts[:,i] = vec(se2vee(b1Tb*bTl))[1:2]
  end
  return kde!( pts, "lcv");
end

function discardOldFeatures!(trkrs::Dict{Int64,Feature})
  for ft in trkrs
    if ft[2].lastzage > 125
      println("deleting $(ft[1])")
      delete!(trkrs, ft[1])
    end
  end
  nothing
end

function propAllTrackers!(trkrs::Dict{Int64,Feature}, bDxb1::Array{Float64,1}, s::Array{Float64,1})
  discardOldFeatures!(trkrs)
  rr = RemoteRef[]
  trkID = Int[]
  for t in trkrs
    # trkrs[t[1]].bel = propagate(t[2].bel, bDxb1, s)
    push!(rr, remotecall(upplocal(), propagate, t[2].bel, bDxb1, s))
    push!(trkID, t[1])

    trkrs[t[1]].age += 1
    trkrs[t[1]].lastzage += 1
  end

  @sync begin
    for i in 1:length(trkID)
      @async trkrs[trkID[i]].bel = fetch(rr[i])
    end
  end

  nothing
end

function p2c(z::Array{Float64,1})
  Rt = R(z[2])
  return Rt*[z[1];0.0], Rt
end

function c2p(x::Array{Float64,1})
  b = atan2(x[2],x[1])
  r = norm(x)
  return r, b
end

function pol2cart(z::Array{Float64,1}, s::Array{Float64,1})
  u, Rt = p2c(z)
  Pp2 = diagm(s.^2)
  P = abs(sqrtm(Rt*Pp2*(Rt')))
  return u, P
end

function cart2pol(z::Array{Float64,1}, s::Array{Float64,1})
  r, b = c2p(z)
  Rt = R(b)
  Pp2 = diagm(s.^2)
  P = abs(sqrtm(Rt'*Pp2*(Rt)))
  return [b;r], P
end

function p2cPtsKDE(z::Array{Float64,1}, s::Array{Float64,1}; N::Int=50)
  u, P = pol2cart(z, s)
  zPts = rand(MvNormal(u,P),N)
  return kde!(zPts, "lcv")
end

featid = map(Int64,0)

function addNewFeatTrk!(trkrs::Dict{Int64,Feature}, z::Array{Float64,1}, s::Array{Float64,1})
  global featid
  # len = length(trkrs)
  featid+=1
  # if length(z) < 3
  #   error("addNewFeat -- z not big enough")
  # end
  trkrs[featid] = Feature(featid, 0, 0, z, p2cPtsKDE(z, s)) #len+1
  return featid # len+1
end

function initTrackersFrom(bearan::Array{Float64,2})
  global featid
  tr = Dict{Int64, Feature}()
  featid = 0
  for i in 1:size(bearan,2)
    addNewFeatTrk!(tr, bearan[:,i], [0.5;0.03])
  end
  return tr
end

function insertLkhdEvals!(lkhds::Array{Float64,2}, currFeats::Dict{Int64,Feature}, sightCart::Array{Float64,2})
  #@show size(lkhds), length(currFeats)
  lkpidx = Array{Int,1}(length(currFeats))
  idx = 0
  for f in currFeats
      yV = evaluateDualTree(f[2].bel, sightCart[1:2,:])
      idx+=1
      lkpidx[idx]=f[1]
      lkhds[:,idx] = yV
  end
  return lkpidx
end

function evalAllLikelihoods(currFeats::Dict{Int64, Feature}, sightFeats::Array{Float64,2})
  numfe = length(currFeats)
  numz = size(sightFeats,2)
  lkhds = zeros(numz, numfe)

  sightCart = zeros(2,numz)
  for i in 1:numz
    u,R = p2c(sightFeats[1:2,i])
    sightCart[:,i] = u
  end

  lkpidx = insertLkhdEvals!(lkhds, currFeats, sightCart)

  return lkhds, lkpidx
end

# agnostic to lkindex permutations
function findNewFeats(lkhds::Array{Float64,2}; thr=0.00001)
  numz = size(lkhds, 1)
  low = vec(maximum(lkhds,2) .< thr)
  return (1:numz)[low]
end

function addNewFeats!(trckr::Dict{Int64,Feature}, lkhds::Array{Float64,2}, lkpidx::Array{Int,1},
                      allmeas::Array{Float64,2}, nidx::Array{Int,1})

  newfeIds = Int[]
  # newfeats = Dict{Int64, Feature}()
  newmeas = allmeas[:,nidx]
  for i in 1:size(newmeas,2)
      id = addNewFeatTrk!(trckr, vec(newmeas[:,i]), [0.5;0.03])
      # newfeats[id] = trckr[id]
      push!(newfeIds, id)
  end
  # increase size of likelihoods
  # r,c = size(lkhds)
  # nlkhds = zeros(r,c+length(newfeIds))
  # nlkhds[:,1:c] = lkhds
  # insertLkhdEvals!(nlkhds, newfeats, allmeas)

  # TODO -- lazy re-evaluation of all likelihoods -- big speedup possible here
  # return evalAllLikelihoods(trckr, allmeas)
  nothing
end

# agnostic to feature idex permutations
function divMaxAcross(lk::Array{Float64,2})
  rlk = round(lk,5)
  m1 = maximum(rlk, 1)
  m1[m1.==0.0] = 1.0
  return rlk./m1
end

# agnostic to featidxpermutations
function divMaxAlong(lk::Array{Float64,2})
  rlk = round(lk,5)
  m2 = maximum(rlk, 2)
  m2[m2.==0.0] = 1.0
  return rlk./m2
end

function hardMatches!(ha::Dict{Int64,Array{Float64,1}}, dmdm::Array{Float64,2}, lkpidx::Array{Int,1}, allmeas::Array{Float64,2})
  # check second maximum is much smaller
  dmdm[dmdm.==2.0]=-1.0
  m1b = (maximum(dmdm,1) .< 0.1)

  for i in 1:size(dmdm,2)
    fid = lkpidx[i]
    if m1b[i] #fid
      zid = findfirst(vec(dmdm[:,i]),-1) # fid
      if zid > 0
        ha[fid] = vec(allmeas[:,zid])
      end
    end
  end
  nothing
end

function findMatches(lk::Array{Float64,2}, lkpidx::Array{Int,1}, allmeas::Array{Float64,2})
  dmac = divMaxAlong(lk)
  dmal = divMaxAcross(lk)
  dmdm = dmac+dmal
  hardassoc = Dict{Int64,Array{Float64,1}}()
  hardMatches!(hardassoc, deepcopy(dmdm), lkpidx, allmeas)
  @show length(hardassoc)
  return hardassoc
end

# will also add new features at the end of the tracking pool
# currently only does hard associations, soft multihypothesis tracking work to follow
function assocMeasWFeats!(trkrs::Dict{Int64, Feature}, fez::Array{Float64,2})
  lk, lkpidx = evalAllLikelihoods(trkrs, fez)
  hardassoc = findMatches(lk, lkpidx, fez)
  # Also add new features into the mix here
  nidx = findNewFeats(lk)
  addNewFeats!(trkrs, lk, lkpidx, fez, nidx)
  return hardassoc
end

function update(bhatXl::BallTreeDensity, z::Array{Float64,1}, s::Array{Float64,1}; N::Int=50)

  bXl = p2cPtsKDE(z,s, N=N)
  # take the product between predicted and measured position
  dummy = kde!(rand(2,N),[1.0])
  pGM, = prodAppxMSGibbsS(dummy, [bhatXl, bXl], Union{}, Union{}, 5)

  # error("update -- Not ready")
  return kde!(pGM, "lcv")
end

function updatelin(bhatXl::BallTreeDensity, z::Array{Float64,1}, s::Array{Float64,1}; N::Int=50)

  bXl = resample(kde!((z')',s),N)
  # take the product between predicted and measured position
  dummy = kde!(rand(length(z),N),[1.0])
  pGM, = prodAppxMSGibbsS(dummy, [bhatXl, bXl], Union{}, Union{}, 5)

  # error("update -- Not ready")
  return kde!(pGM, "lcv")
end

function updatelin(bhatXl::BallTreeDensity, bXl::BallTreeDensity; N=50)
  dummy = kde!(rand(bhatXl.bt.dims,N),[1.0])
  pGM, = prodAppxMSGibbsS(dummy, [bhatXl, bXl], Union{}, Union{}, 5)
  # error("update -- Not ready")
  return kde!(pGM, "lcv")
end

function doAsyncUpdate(oldFeat::Feature, meas::Array{Float64,1}, s::Array{Float64,1})
  bel = update(oldFeat.bel, meas, s)
  return Feature(oldFeat.id, oldFeat.age, 0, meas[1:3], bel)
end

function measUpdateTrackers!(trkrs::Dict{Int64,Feature}, assoc::Dict{Int64,Array{Float64,1}}, s::Array{Float64,1})
  # loop through all associated measurements
  # TODO --this can all be done in parallel to reduce computation time
  # dbg = Dict{Int64,Array{BallTreeDensity,1}}()
  rr = RemoteRef[]
  trkID = Int[]
  for meas in assoc
    # bel = deepcopy(trkrs[meas[1]].bel)
    # trkrs[meas[1]].bel = update(bel, meas[2], s)
    # trkrs[meas[1]].lastzage = 0
    # trkrs[meas[1]].lastz = meas[2]

    # trkrs[meas[1]] = doAsyncUpdate(trkrs[meas[1]], meas[2], s)
    push!(rr, remotecall(upplocal(), doAsyncUpdate, trkrs[meas[1]], meas[2], s))
    push!( trkID, meas[1])

    # arr = Array{BallTreeDensity,1}()
    # push!(arr,bel)
    # push!(arr,p2cPtsKDE(meas[2],[1.0;0.02]))
    # push!(arr,deepcopy(trkrs[meas[1]].bel))
    # dbg[meas[1]] = arr
  end

  @sync begin
    for i in 1:length(trkID)
      @async trkrs[trkID[i]] = fetch(rr[i])
    end
  end

  # return dbg
  nothing
end
