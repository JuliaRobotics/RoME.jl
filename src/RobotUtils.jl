
# warn("Must remove optional dev/ISAMRemoteSolve.jl dependency here.")
# include("dev/ISAMRemoteSolve.jl")


immutable RangeAzimuthElevation
  range::Float64
  azimuth::Float64
  elevation::Union{Void,Float64}
end

function convert{T <: CoordinateTransformations.AffineMap}(::Type{T}, x::SE3)
  q = convert(TransformUtils.Quaternion, x.R)
  Translation(x.t...) ∘ LinearMap( Rotations.Quat(q.s, q.v...) )
end
function convert{T <: CoordinateTransformations.AffineMap{Rotations.Quat{Float64}}}(::Type{SE3}, x::T)
  SE3(x.v[1:3], TransformUtils.Quaternion(x.m.w, [x.m.x,x.m.y,x.m.z]) )
end

function convert(::Type{SE3}, t::Tuple{Symbol, Vector{Float64}})
  if t[1]==:XYZqWXYZ
    return SE3(t[2][1:3],TransformUtils.Quaternion(t[2][4],t[2][5:7]))
  else
    error("Unknown conversion type $(t[1])")
  end
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

# should be deprecated or indicated more clearly
lsrBR(a) = [a[2,:];a[1,:]]';

function veePose3(s::SE3)
  TransformUtils.veeEuler(s)
end
function veePose(s::SE3)
  TransformUtils.veeEuler(s)
end


function \(s::SE3, wTr::CTs.Translation)
  bTr = s.R.R'*(wTr.v-s.t)
  Dtr = bTr
  range = norm(Dtr)
  azi = atan2(Dtr[2], Dtr[1])
  elev = atan2(Dtr[3], Dtr[1])
  RangeAzimuthElevation(range, azi, elev)
end



"""
    getRangeKDEMax2D(fgl::FactorGraph, vsym1::Symbol, vsym2::Symbol)

Calculate the cartesian distance between two vertices in the graph using their symbol name, and by maximum belief point.
"""
function getRangeKDEMax2D(fgl::FactorGraph, vsym1::Symbol, vsym2::Symbol)
  x1 = getKDEMax(getVertKDE(fgl, vsym1))
  x2 = getKDEMax(getVertKDE(fgl, vsym2))
  norm(x1[1:2]-x2[1:2])
end

function measureMeanDist(fg::FactorGraph, a::T, b::T) where {T <: AbstractString}
    #bearrang!(residual::Array{Float64,1}, Z::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
    res = zeros(2)
    A = getVal(fg,a)
    B = getVal(fg,b)
    Ax = Base.mean(vec(A[1,:]))
    Ay = Base.mean(vec(A[2,:]))
    Bx = Base.mean(vec(B[1,:]))
    By = Base.mean(vec(B[2,:]))
    dx = Bx - Ax
    dy = By - Ay
    b = atan2(dy,dx)
    r = sqrt(dx^2 + dy^2)
    return r, b
end

function predictBodyBR(fg::FactorGraph, a::T, b::T) where {T <: AbstractString}
  res = zeros(2)
  A = getVal(fg,a)
  B = getVal(fg,b)
  Ax = Base.mean(vec(A[1,:]))
  Ay = Base.mean(vec(A[2,:]))
  Ath = Base.mean(vec(A[3,:]))
  Bx = Base.mean(vec(B[1,:]))
  By = Base.mean(vec(B[2,:]))
  wL = SE2([Bx;By;0.0])
  wBb = SE2([Ax;Ay;Ath])
  bL = se2vee((wBb \ eye(3)) * wL)
  dx = bL[1] - 0.0
  dy = bL[2] - 0.0
  b = (atan2(dy,dx))
  r = sqrt(dx^2 + dy^2)
  return b, r
end

function getNextLbl(fgl::FactorGraph, chr)
  # TODO convert this to use a double lookup
  warn("getNextLbl(::FactorGraph..) to be deprecated, use getlastpose/landm(::SLAMWrapper..) instead.")
  max = -1
  maxid = -1
  for vid in fgl.IDs
  # for v in fgl.v #fgl.g.vertices # fgl.v
      v = (vid[2], fgl.g.vertices[vid[2]])
      if v[2].attributes["label"][1] == chr
        val = parse(Int,v[2].attributes["label"][2:end])
        if max < val
          max = val
          maxid = v[1]
        end
      end
  end
  if maxid != -1
    v = getVert(fgl,maxid)
    X = getVal(v)
    return v, X, Symbol(string(chr,max+1))
  else
    return nothing, nothing, Symbol(string(chr,max+1)) # Union{}
  end
end

function getLastPose(fgl::FactorGraph)
  return getNextLbl(fgl, 'x')
end
getLastPose2D(fgl::FactorGraph) = getLastPose(fgl)

function getlastpose(slam::SLAMWrapper)
  error("getlastpose -- Not implemented yet")
end


function getLastLandm2D(fgl::FactorGraph)
  return getNextLbl(fgl, 'l')
end

function odomKDE(p1,dx,cov)
  X = getPoints(p1)
  sig = diag(cov)
  RES = zeros(size(X))
  # increases the number of particles based on the number of modes in the measurement Z
  for i in 1:size(X,2)
      ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
      RES[:,i] = addPose2Pose2(X[:,i], dx + ent)
  end
  return kde!(RES, "lcv")
end

"""
    addOdoFG!(fg, name, DX, cov, , N=100, ready=1)

Create a new variable node and insert odometry constraint factor between
which will automatically increment latest pose symbol x<k+1> for new node new node and
constraint factor are returned as a tuple.

"""
function addOdoFG!(
        fg::FactorGraph,
        n::Symbol,
        DX::Array{Float64,1},
        cov::Array{Float64,2};
        N::Int=0,
        ready::Int=1,
        labels::Vector{<:AbstractString}=String[]  )
    #
    prev, X, nextn = getLastPose2D(fg)
    r,c = size(X)
    if N==0
      N = c
    end
    sig = diag(cov)
    XnextInit = zeros(r,c)
    # increases the number of particles based on the number of modes in the measurement Z
    for i in 1:c
        ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
        XnextInit[:,i] = addPose2Pose2(X[:,i], DX + ent)
    end

    v = addNode!(fg, n, XnextInit, cov, N=N, ready=ready, labels=labels)
    pp = Pose2Pose2(vectoarr2(DX), cov, [1.0]) #[prev;v],
    f = addFactor!(fg, [prev;v], pp, ready=ready)
    infor = inv(cov^2)
    # addOdoRemote(prev.index,v.index,DX,infor) # this is for remote factor graph ref parametric solution -- skipped internally by global flag variable
    return v, f
end

function addOdoFG!(
        fg::FactorGraph,
        n::T,
        Z::Pose3Pose3;
        N::Int=0,
        ready::Int=1,
        labels::Vector{T}=String[]  ) where {T <: AbstractString}
  #
  DX=Z.μ
  cov=Z.Σ.mat

  addOdoFG(fg, n, DX, cov, N=N, ready=ready, labels=labels)
end

"""
    addOdoFG!(fg, odo, N=100, ready=1)

Create a new variable node and insert odometry constraint factor between
which will automatically increment latest pose symbol x<k+1> for new node new node and
constraint factor are returned as a tuple.

"""
function addOdoFG!(
        fgl::FactorGraph,
        odo::PP;
        N::Int=0,
        ready::Int=1,
        labels::Vector{T}=String[] ) where {PP <: RoME.BetweenPoses, T <: AbstractString}
    #
    vprev, X, nextn = getLastPose(fgl)
    if N==0
      N = size(X,2)
    end
    vnext = addNode!(fgl, nextn, X⊕odo, ones(1,1), N=N, ready=ready, labels=labels)
    fact = addFactor!(fgl, [vprev;vnext], odo)
    return vnext, fact
end


function initfg(;sessionname="NA")
  fgl = emptyFactorGraph()
  fgl.sessionname=sessionname
  # registerCallback!(fgl, RoME.getSample) # RoME.evalPotention
  return fgl
end

function initFactorGraph!(fg::FactorGraph;
      P0::Union{Array{Float64,2},Void}=nothing,
      init::Union{Vector{Float64},Void}=nothing,
      N::Int=100,
      lbl::Symbol=:x0,
      ready::Int=1,
      firstPoseType=Pose2,
      labels::Vector{<:AbstractString}=String[])
  #
  nodesymbols = Symbol[]
  if firstPoseType == Pose2
      init = init!=nothing ? init : zeros(3)
      P0 = P0!=nothing ? P0 : diagm([0.03;0.03;0.001])
      init = vectoarr2(init)
      addNode!(fg,lbl,Pose2,N=N,autoinit=true,ready=ready,labels=String["VARIABLE"; labels] )
      push!(nodesymbols, lbl)
      # v1 = addNode!(fg, lbl, init, P0, N=N, ready=ready, labels=labels)
      fctVert = addFactor!(fg, [lbl;], PriorPose2(init, P0,  [1.0]), ready=ready, labels=String["FACTOR"; labels]) #[v1],
      push!(nodesymbols, fctVert.label)
  end
  if firstPoseType == Pose3
      init = init!=nothing ? init : zeros(6)
      P0 = P0!=nothing ? P0 : diagm([0.03;0.03;0.03;0.001;0.001;0.001])
      addNode!(fg,lbl,Pose2,N=N,autoinit=true,ready=ready,labels=String["VARIABLE"; labels] )
      push!(nodesymbols, lbl)
      # v1 = addNode!(fg, lbl, init, P0, N=N, ready=ready, labels=labels)
      fctVert = addFactor!(fg, [lbl;], PriorPose3(MvNormal(init, P0)), ready=ready, labels=String["FACTOR"; labels]) #[v1],
      push!(nodesymbols, fctVert.label)
  end
  warn("initFactorGraph! -- returns variable and factor symbols ")
  return nodesymbols
end


function newLandm!(fg::FactorGraph, lm::T, wPos::Array{Float64,2}, sig::Array{Float64,2};
                  N::Int=100, ready::Int=1, labels::Vector{T}=String[]) where {T <: AbstractString}

    # TODO -- need to confirm this function is updating the correct memory location. v should be pointing into graph
    vert=addNode!(fg, Symbol(lm), wPos, sig, N=N, ready=ready, labels=labels)

    vert.attributes["age"] = 0
    vert.attributes["maxage"] = 0
    vert.attributes["numposes"] = 0
    updateFullVert!(fg, vert)

    println("newLandm! -- added $(lm)")
    return vert
end

function updateLandmAge(vlm::Graphs.ExVertex, pose::T) where {T <: AbstractString}
  error("still working here")
end

function addBRFG!(fg::FactorGraph,
      pose::T,
      lm::T,
      br::Array{Float64,1},
      cov::Array{Float64,2};
      ready::Int=1  ) where {T <: AbstractString}
  #
  vps = getVert(fg,pose)
  vlm = getVert(fg,lm)
  testlbl = vps.label*vlm.label
  for nei in getOutNeighbors(fg, vlm)
    if nei.label == testlbl
      # TODO -- makes function call brittle
      warn("We already have $(testlbl), skipping this constraint")
      return nothing
    end
  end
  @show keys(vlm.attributes)
  np = vlm.attributes["numposes"]
  la = vlm.attributes["age"]
  nage = parse(Int,pose[2:end])
  vlm.attributes["numposes"] = np+1
  vlm.attributes["age"] = ((la*np)+nage)/(np+1)
  vlm.attributes["maxage"] = nage
  updateFullVert!(fg, vlm)

  # @show br, cov
  pbr = Pose2DPoint2DBearingRange{Normal, Normal}(Normal(br[1], cov[1,1]), Normal(br[2],  cov[2,2]))
  # pbr = Pose2DPoint2DBearingRange((br')',  cov,  [1.0])
  f = addFactor!(fg, [vps;vlm], pbr, ready=ready ) #[vps;vlm],

  # only used for max likelihood unimodal tests.
  u, P = pol2cart(br[[2;1]], diag(cov))
  infor = inv(P^2)
  # addLandmMeasRemote(vps.index,vlm.index,u,infor) # for iSAM1 remote solution as reference
  return f
end

function addMMBRFG!(fg::FactorGraph, pose::T,
                  lm::Array{T,1}, br::Array{Float64,1},
                  cov::Array{Float64,2}; w=[0.5;0.5], ready::Int=1) where {T <: AbstractString}
    #
    vps = getVert(fg,pose)
    vlm1 = getVert(fg,lm[1])
    vlm2 = getVert(fg,lm[2])

    pbr = Pose2DPoint2DBearingRange(vectoarr2(br),  cov,  w) #[vps;vlm1;vlm2],
    f = addFactor!(fg, [vps;vlm1;vlm2], pbr, ready=ready )
    return f
end


function projNewLandmPoints(vps::Graphs.ExVertex, br::Array{Float64,1}, cov::Array{Float64,2})
    # TODO -- convert to use Distributions and common projection function
    Xps = getVal(vps)
    lmPts = zeros(2,size(Xps,2))
    for i in 1:size(Xps,2)
        ent = [cov[1,1]*randn(); cov[2,2]*randn()]
        init = vec(Xps[1:2,i])+randn(2)
        lmPts[:,i] = solveLandm(br + ent, vec(Xps[:,i]), init)
    end
    return lmPts
end

function projNewLandm!(fg::FactorGraph, pose::T, lm::T, br::Array{Float64,1}, cov::Array{Float64,2};
                        addfactor=true, N::Int=100, ready::Int=1,
                        labels::Vector{T}=String[]) where {T <: AbstractString}
    #
    vps = getVert(fg,pose)

    lmPts = projNewLandmPoints(vps, br, cov)
    vlm = newLandm!(fg, lm, lmPts, cov, N=N, ready=ready, labels=labels) # cov should not be required here
    if addfactor
      fbr = addBRFG!(fg, pose, lm, br, cov, ready=ready)
      return vlm, fbr
    end
    return vlm
end

function calcIntersectVols(fgl::FactorGraph, predLm::BallTreeDensity;
                          currage=0, maxdeltaage=Inf)
    # all landmarks of interest
    xx,ll = ls(fgl)
    # output result
    rr = Dict{String, RemoteRef}()
    fetchlist = String[]
    iv = Dict{String, Float64}()
    for l in ll
      pvlm = getVert(fgl,l)
      # TODO -- can be improved via query in DB case
      if currage - pvlm.attributes["maxage"] < maxdeltaage
        p = getVertKDE(fgl, l)
        rr[l] = remotecall(uppA(), intersIntgAppxIS, p,predLm)
        push!(fetchlist, l)
      else
        println("calcIntersectVols -- ignoring $(l) because maxdeltaage exceeded")
        iv[l] = 0
      end
    end
    max = 0
    maxl = String("")
    for l in fetchlist #ll
      # p = getVertKDE(fgl, l)
      # tv = intersIntgAppxIS(p,predLm)
      # iv[l] = tv
      tv = fetch(rr[l])
      iv[l] = tv
      if max < tv  max = tv; maxl = l; end
    end
    return iv, maxl
end

function maxIvWithoutID(ivs::Dict{String, Float64}, l::T) where {T <: AbstractString}
  max = 0
  maxl = String("")
  for i in ivs
    if max < i[2] && i[1] != l;  max = i[2]; maxl = i[1]; end
  end
  return maxl
end

# binary tests to distinguish how to automatically add a landmark to the existing factor graph
function doAutoEvalTests(fgl::FactorGraph, ivs::Dict{T, Float64}, maxl::T, lmid::Int, lmindx::Int) where {T <: AbstractString}
  maxAnyval = maxl != String("") ? ivs[maxl] : 0.0
  # maxid = fgl.IDs[maxl]
  lmidSugg = lmid != -1 # a landmark ID has been suggested
  maxAnyExists = maxAnyval > 0.03 # there is notable intersection with previous landm
  lmIDExists = haskey(fgl.v, lmid) # suggested lmid already in fgl
  newlmindx = lmindx
  if lmIDExists
    lmSuggLbl = String(getVert(fgl,lmid).label) # TODO -- wasteful
  else
    newlmindx = lmindx + 1
    lmSuggLbl = String(string('l',newlmindx))
  end
  maxl2 = lmIDExists ? maxIvWithoutID(ivs, lmSuggLbl) : String("")
  maxl2Exists = lmIDExists ? (maxl2 != "" ? ivs[maxl2] > 0.03 : false) : false # there is notable intersection with previous landm
  intgLmIDExists = lmIDExists ? ivs[lmSuggLbl] > 0.03 : false

  return lmidSugg, maxAnyExists, maxl2Exists, maxl2, lmIDExists, intgLmIDExists, lmSuggLbl, newlmindx
end


function evalAutoCases!(fgl::FactorGraph, lmid::Int, ivs::Dict{T, Float64}, maxl::T,
                        pose::T, lmPts::Array{Float64,2}, br::Array{Float64,1}, cov::Array{Float64,2}, lmindx::Int;
                        N::Int=100, ready::Int=1 ) where {T <: AbstractString}
  lmidSugg, maxAnyExists, maxl2Exists, maxl2, lmIDExists, intgLmIDExists, lmSuggLbl, newlmindx = doAutoEvalTests(fgl,ivs,maxl,lmid, lmindx)

  println("evalAutoCases -- found=$(lmidSugg), $(maxAnyExists), $(maxl2Exists), $(lmIDExists), $(intgLmIDExists)")

  vlm = Union{}; fbr = Union{};
  if (!lmidSugg && !maxAnyExists)
    #new landmark and UniBR constraint
    v,L,lm = getLastLandm2D(fgl)
    vlm = newLandm!(fgl, lm, lmPts, cov, N=N,ready=ready)
    fbr = addBRFG!(fgl, pose, lm, br, cov, ready=ready)
  elseif !lmidSugg && maxAnyExists
    # add UniBR to best match maxl
    vlm = getVert(fgl,maxl)
    fbr = addBRFG!(fgl, pose, maxl, br, cov, ready=ready)
  elseif lmidSugg && !maxl2Exists && !lmIDExists
    #add new landmark and add UniBR to suggested lmid
    vlm = newLandm!(fgl, lmSuggLbl, lmPts, cov, N=N, ready=ready)
    fbr = addBRFG!(fgl, pose, lmSuggLbl, br, cov, ready=ready)
  elseif lmidSugg && !maxl2Exists && lmIDExists && intgLmIDExists
    # doesn't self intesect with existing lmid, add UniBR to lmid
    vlm = getVert(fgl, lmid)
    fbr = addBRFG!(fgl, pose, lmSuggLbl, br, cov, ready=ready)
  elseif lmidSugg && maxl2Exists && !lmIDExists
    # add new landmark and add MMBR to both maxl and lmid
    vlm = newLandm!(fgl, lmSuggLbl, lmPts, cov, N=N, ready=ready)
    addMMBRFG!(fgl, pose, [maxl2;lmSuggLbl], br, cov, ready=ready)
  elseif lmidSugg && maxl2Exists && lmIDExists && intgLmIDExists
    # obvious case, add MMBR to both maxl and lmid. Double intersect might be the same thing
    println("evalAutoCases! -- obvious case is happening")
    addMMBRFG!(fgl, pose, [maxl2;lmSuggLbl], br, cov, ready=ready)
    vlm = getVert(fgl,lmSuggLbl)
  elseif lmidSugg && maxl2Exists && lmIDExists && !intgLmIDExists
    # odd case, does not intersect with suggestion, but does with some previous landm
    # add MMBR
    warn("evalAutoCases! -- no self intersect with suggested $(lmSuggLbl) detected")
    addMMBRFG!(fgl, pose, [maxl;lmSuggLbl], br, cov, ready=ready)
    vlm = getVert(fgl,lmSuggLbl)
  elseif lmidSugg && !maxl2Exists && lmIDExists && !intgLmIDExists
  #   # landm exists but no intersection with existing or suggested lmid
  #   # may suggest some error
    warn("evalAutoCases! -- no intersect with suggested $(lmSuggLbl) or map detected, adding  new landmark MM constraint incase")
    v,L,lm = getLastLandm2D(fgl)
    vlm = newLandm!(fgl, lm, lmPts, cov, N=N, ready=ready)
    addMMBRFG!(fgl, pose, [lm; lmSuggLbl], br, cov, ready=ready)
  else
    error("evalAutoCases! -- unknown case encountered, can reduce to this error to a warning and ignore user request")
  end

  return vlm, fbr, newlmindx
end

function addAutoLandmBR!(fgl::FactorGraph, pose::T, lmid::Int, br::Array{Float64,1}, cov::Array{Float64,2}, lmindx::Int;
                      N::Int=100, ready::Int=1 ) where {T <: AbstractString}
    vps = getVert(fgl, pose)
    lmPts = projNewLandmPoints(vps, br, cov)
    lmkde = kde!(lmPts)
    currage = parse(Int, pose[2:end])
    ivs, maxl = calcIntersectVols(fgl, lmkde, currage=currage,maxdeltaage=10)

    # There are 8 cases of interest
    vlm, fbr, newlmindx = evalAutoCases!(fgl, lmid, ivs, maxl,pose,lmPts, br,cov,lmindx,N=N,ready=ready)

    return vlm, fbr, newlmindx
end

function malahanobisBR(measA, preA, cov::Array{Float64,2})
    # measure landmark with noise
    res = measA - preA
    mala2 = Union{}
    #Malahanobis distance
    if false
      lambda = cov \ eye(2)
      mala2 = res' * lambda * res
    else
      mala2 = res' * (cov \ res)
    end

    mala = sqrt(mala2)
    return mala
end





# ------------------------------------
# Transfered from IncrementalInference



function get2DSamples(fg::FactorGraph,
      sym;
      from::Int=0, to::Int=999999999,
      minnei::Int=0,
      api::DataLayerAPI=IncrementalInference.localapi  )
  #
  X = Array{Float64,1}()
  Y = Array{Float64,1}()

  # if sym = 'l', ignore single measurement landmarks

  for id in fg.IDs
    vertlbl = string(id[1])
    if vertlbl[1] == sym
      val = parse(Int,vertlbl[2:end])
      if from <= val && val <= to
        if length( getOutNeighbors(fg, id[2], api=api, needdata=true ) ) >= minnei
          # if length(out_neighbors(fg.v[id[2]],fg.g)) >= minnei
          X=[X; vec(getVal(fg,id[2], api=api)[1,:]) ]
          Y=[Y; vec(getVal(fg,id[2], api=api)[2,:]) ]
          # X=[X;vec(fg.v[id[2]].attributes["val"][1,:])]
          # Y=[Y;vec(fg.v[id[2]].attributes["val"][2,:])]
        end
      end
    end
  end
  return X,Y
end

function getAll2D(fg, sym; minnei::Int=0, api::DataLayerAPI=IncrementalInference.localapi)
  warn("getAll2D deprecated, use get2DSamples instead")
  return get2DSamples(fg, sym, minnei=minnei, api=api)
end

function get2DSampleMeans(fg::FactorGraph,
      sym;
      from::Int=0, to::Int=9999999999,
      minnei::Int=0,
      api::DataLayerAPI=IncrementalInference.localapi  )
  #
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]

  # if sym = 'l', ignore single measurement landmarks
  allIDs = Array{Int,1}()
  for id in fg.IDs
    vertlbl = string(id[1])
    if vertlbl[1] == sym
      val = parse(Int,vertlbl[2:end])
      if from <= val && val <= to
        if length( getOutNeighbors(fg, id[2], api=api, needdata=true) ) >= minnei
          push!(allIDs, val)
        end
      end
    end
  end
  allIDs = sort(allIDs)

  for id in allIDs
    X=[X;Base.mean( vec( getVal(fg,string(sym,id), api=api)[1,:] ) )]
    Y=[Y;Base.mean( vec( getVal(fg, string(sym,id), api=api)[2,:] ) )]
    if sym == 'x'
      Th=[Th;Base.mean( vec( getVal(fg, string(sym,id), api=api)[3,:] ) )]
    end
    push!(LB, string(sym,id))
  end
  return X,Y,Th,LB
end

#draw landmark positions
function getAll2DMeans(fg, sym, api::DataLayerAPI=IncrementalInference.localapi)
  return get2DSampleMeans(fg, sym, api=api)
end

function getAll2DPoses(fg::FactorGraph, api::DataLayerAPI=IncrementalInference.localapi)
    return getAll2DSamples(fg, 'x', api=api)
end

function get2DPoseSamples(fg::FactorGraph; from::Int=0, to::Int=999999999, api::DataLayerAPI=IncrementalInference.localapi)
  return get2DSamples(fg, 'x'; from=from, to=to, api=api)
end

function get2DPoseMeans(fg::FactorGraph; from::Int=0, to::Int=999999999, api::DataLayerAPI=IncrementalInference.localapi)
  return get2DSampleMeans(fg, 'x', from=from, to=to, api=api)
end


function get2DPoseMax(fgl::FactorGraph;
            from::Int=-99999999999, to::Int=9999999999, api::DataLayerAPI=IncrementalInference.localapi )
  xLB,ll = ls(fgl) # TODO add: from, to, special option 'x'
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]
  for slbl in xLB
    lbl = string(slbl)
    if from <= parse(Int,lbl[2:end]) <=to
      mv = getKDEMax(getVertKDE(fgl,slbl))
      push!(X,mv[1])
      push!(Y,mv[2])
      push!(Th,mv[3])
      push!(LB, lbl)
    end
  end
  return X, Y, Th, LB
end

function getAll2DLandmarks(fg::FactorGraph, minnei::Int=0, api::DataLayerAPI=IncrementalInference.localapi)
    return getAll2DSamples(fg, 'l', minnei=minnei, api=api)
end

function get2DLandmSamples(fg::FactorGraph; from::Int=0, to::Int=999999999, minnei::Int=0, api::DataLayerAPI=IncrementalInference.localapi)
  return get2DSamples(fg, 'l', from=from, to=to, minnei=minnei, api=api)
end

function get2DLandmMeans(fg::FactorGraph; from::Int=0, to::Int=999999999, minnei::Int=0, api::DataLayerAPI=IncrementalInference.localapi)
  return get2DSampleMeans(fg, 'l', from=from, to=to, minnei=minnei, api=api)
end

function removeKeysFromArr(fgl::FactorGraph, torm::Array{Int,1}, lbl::Array{String,1})
  retlbs = String[]
  for i in 1:length(lbl)
    id = parse(Int,lbl[i][2:end])
    if findfirst(torm,id) == 0
      push!(retlbs, lbl[i])
    else
      println("removeKeysFromArr -- skipping $(lbl[i]), id=$(id)")
    end
  end
  return retlbs
end

function get2DLandmMax(fgl::FactorGraph;
                from::Int=-99999999999,
                to::Int=9999999999,
                showmm=false, MM=Union{},
                api::DataLayerAPI=IncrementalInference.localapi  )
  #
  xLB,lLB = ls(fgl) # TODO add: from, to, special option 'x'
  if !showmm lLB = removeKeysFromArr(fgl, collect(keys(MM)), lLB); end
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]
  for lb in lLB
    @show lb
    @show lbl = string(lb)
    if from <= parse(Int,lbl[2:end]) <=to
      mv = getKDEMax(getVertKDE(fgl,lb, api=api))
      push!(X,mv[1])
      push!(Y,mv[2])
      push!(LB, lbl)
    end
  end
  return X, Y, Th, LB
end



# convenience function to add DIDSON sonar constraints to graph
function addLinearArrayConstraint(fgl::FactorGraph,
      rangebearing::Union{Tuple{Float64, Float64}, Vector{Float64}},
      pose::Symbol,
      landm::Symbol ;
      rangecov::Float64=3e-4,
      bearingcov::Float64=3e-4 )

  #

  cl = LinearRangeBearingElevation((rangebearing[1],rangecov),(rangebearing[2],bearingcov))
  if !haskey(fgl.IDs, landm)
    pts = getVal(fgl, pose) + cl
    N = size(pts,2)
    vl1 = addNode!(fgl, landm,  pts,  N=N)
    println("Automatically added $(landm) to the factor graph")
  end
  addFactor!(fgl, [getVert(fgl, pose); getVert(fgl, landm)], cl)
  nothing
end


function addSoftEqualityPoint2D(fgl::FactorGraph, l1::Symbol, l2::Symbol;
    dist=MvNormal([0.0;0.0],eye(2)), ready::Int=1 )
  pp = Point2DPoint2D(dist)
  addFactor!(fgl, [l1,l2], pp, ready=ready)
end



#
# function +(p1::Rigid6DOF, p2::Rigid6DOF)
#     t = rotate(p1.rot, p2.trl) + p1.trl
#     return Rigid6DOF(p1.rot*p2.rot,t)
# end
#
# function +(p1::PoseSE3, p2::PoseSE3)
#     return PoseSE3(p1.utime+p2.utime,
#                    "+",
#                    p1.mu+p2.mu,
#                    p1.cov + p2.cov,
#                    Dict())
# end
#
# function makePose3(q::Quaternion=Quaternion(1.,zeros(3)), t::Array{Float64,1}=[0.,0,0]; ut=0, name="pose", cov=0.1*eye(6))
#     return PoseSE3(ut,
#                    name,
#                    Rigid6DOF( q, t ),
#                    cov,
#                    Dict())
# end

# function makeDynPose3(q::Quaternion=Quaternion(1.,zeros(3)), t::Array{Float64,1}=[0.,0,0]; ut=0, name="pose", cov=0.1*eye(18))
#     return DynPose3(ut,
#                     name,
#                     Dynamic6DOF( q, t, [0.,0,0], [0.,0,0] ),
#                     IMUComp([0.,0,0],[0.,0,0]),
#                     cov,
#                     Dict())
# end
#
# function wTo(p::PoseSE3)
#     T = eye(4)
#     T[1:3,1:3] =  convert(SO3, p.mu.rot).R
#     T[1:3,4] = p.mu.trl
#     return T
# end
#
# function oTw(p::PoseSE3)
#     return wTo(p) \ eye(4)
# end
