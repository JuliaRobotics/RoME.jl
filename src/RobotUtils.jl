
include("dev/ISAMRemoteSolve.jl")

function measureMeanDist(fg::FactorGraph, a::ASCIIString, b::ASCIIString)
    #bearrang!(residual::Array{Float64,1}, Z::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
    res = zeros(2)
    A = getVal(fg.v[fg.IDs[a]])
    B = getVal(fg.v[fg.IDs[b]])
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

function predictBodyBR(fg::FactorGraph, a::ASCIIString, b::ASCIIString)
  res = zeros(2)
  A = getVal(fg.v[fg.IDs[a]])
  B = getVal(fg.v[fg.IDs[b]])
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

# function getLastPose2D(fg::FactorGraph)
#   max = 0
#   for id in fg.IDs
#       if id[1][1] == 'x'
#         val = parse(Int,id[1][2:end])
#         if max < val
#           max = val
#         end
#       end
#   end
#   v = fg.v[fg.IDs[string('x',max)]]
#   X = getVal(v)
#   return v, X, string('x',max+1)
# end

function getNextLbl(fgl::FactorGraph, sym)
  max = 0
  maxid = -1
  for v in fgl.v
      if v[2].attributes["label"][1] == sym
        val = parse(Int,v[2].attributes["label"][2:end])
        if max < val
          max = val
          maxid = v[1]
        end
      end
  end
  if maxid != -1
    v = fgl.v[maxid]
    X = getVal(v)
    return v, X, string(sym,max+1)
  else
    return Union{}, Union{}, string(sym,max+1)
  end
end

function getLastPose2D(fgl::FactorGraph)
  return getNextLbl(fgl, 'x')
  # max = 0
  # maxid = 0
  # for v in fg.v
  #     if v[2].attributes["label"][1] == 'x'
  #       val = parse(Int,v[2].attributes["label"][2:end])
  #       if max < val
  #         max = val
  #         maxid = v[1]
  #       end
  #     end
  # end
  # v = fg.v[maxid]
  # X = getVal(v)
  # return v, X, string('x',max+1)
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


function addOdoFG!(fg::FactorGraph, n::ASCIIString, DX::Array{Float64,1}, cov::Array{Float64,2};
                  N::Int=100)
    prev, X, nextn = getLastPose2D(fg)
    sig = diag(cov)
    RES = zeros(size(X))
    # increases the number of particles based on the number of modes in the measurement Z
    for i in 1:size(X,2)
        ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
        RES[:,i] = addPose2Pose2(X[:,i], DX + ent)
    end

    v = addNode!(fg, n, RES, cov, N=N)
    pp = Pose2Pose2([prev;v], (DX')', cov, [1.0])
    f = addFactor!(fg, pp)
    infor = inv(cov^2)
    addOdoRemote(prev.index,v.index,DX,infor) # this is for remote factor graph ref parametric solution -- skipped internally by global flag variable
    return v, f
end

function initfg()
  return emptyFactorGraph()
end


function newLandm!(fg::FactorGraph, lm::ASCIIString, wPos::Array{Float64,2}, sig::Array{Float64,2};
                  N::Int=100)
    v=addNode!(fg, lm, wPos, sig, N=N)
    fg.v[v.index].attributes["age"] = 0
    fg.v[v.index].attributes["maxage"] = 0
    fg.v[v.index].attributes["numposes"] = 0
    println("newLandm! -- added $(lm)")
    return v
end

function updateLandmAge(vlm::Graphs.ExVertex, pose::AbstractString)

end

function addBRFG!(fg::FactorGraph, pose::ASCIIString,
                  lm::ASCIIString, br::Array{Float64,1}, cov::Array{Float64,2})
    vps = fg.v[fg.IDs[pose]]
    vlm = fg.v[fg.IDs[lm]]
    # println("addBRFG! -- looking at $(lm)")
    np = vlm.attributes["numposes"]
    la = vlm.attributes["age"]
    nage = parse(Int,pose[2:end])
    vlm.attributes["numposes"] = np+1
    vlm.attributes["age"] = ((la*np)+nage)/(np+1)
    vlm.attributes["maxage"] = nage
    f = addFactor!(fg, Pose2DPoint2DBearingRange([vps;vlm],(br')',  cov,  [1.0]) )

    # only used for max likelihood unimodal tests.
    u, P = pol2cart(br[[2;1]], diag(cov))
    infor = inv(P^2)
    addLandmMeasRemote(vps.index,vlm.index,u,infor) # for iSAM1 remote solution as reference
    return f
end

function addMMBRFG!(fg::FactorGraph, pose::ASCIIString,
                  lm::Array{ASCIIString,1}, br::Array{Float64,1},
                  cov::Array{Float64,2}; w=[0.5;0.5])
    vps = fg.v[fg.IDs[pose]]
    vlm1 = fg.v[fg.IDs[lm[1]]]
    vlm2 = fg.v[fg.IDs[lm[2]]]

    pbr = Pose2DPoint2DBearingRange([vps;vlm1;vlm2],(br')',  cov,  w)
    f = addFactor!(fg, pbr )
    return f
end

function projNewLandmPoints(vps::Graphs.ExVertex, br::Array{Float64,1}, cov::Array{Float64,2})
    Xps = getVal(vps)
    lmPts = zeros(2,size(Xps,2))
    for i in 1:size(Xps,2)
        ent = [cov[1,1]*randn(); cov[2,2]*randn()]
        init = vec(Xps[1:2,i])+randn(2)
        lmPts[:,i] = solveLandm(br + ent, vec(Xps[:,i]), init)
    end
    return lmPts
end

function projNewLandm!(fg::FactorGraph, pose::ASCIIString, lm::ASCIIString, br::Array{Float64,1}, cov::Array{Float64,2};
                        addfactor=true, N::Int=100)
    vps = fg.v[fg.IDs[pose]]
    # Xps = getVal(vps)
    # lmPts = zeros(2,size(Xps,2))
    # for i in 1:size(Xps,2)
    #     ent = [cov[1,1]*randn(); cov[2,2]*randn()]
    #     init = vec(Xps[1:2,i])+randn(2)
    #     lmPts[:,i] = solveLandm(br + ent, vec(Xps[:,i]), init)
    # end
    lmPts = projNewLandmPoints(vps, br, cov)
    vlm = newLandm!(fg, lm, lmPts, cov, N=N) # cov should not be required here
    if addfactor
      fbr = addBRFG!(fg, pose, lm, br, cov)
      return vlm, fbr
    end
    return vlm
end

function calcIntersectVols(fgl::FactorGraph, predLm::BallTreeDensity;
                          currage=0, maxdeltaage=Inf)
    # all landmarks of interest
    xx,ll = ls(fgl)
    # output result
    rr = Dict{ASCIIString, RemoteRef}()
    fetchlist = ASCIIString[]
    iv = Dict{ASCIIString, Float64}()
    for l in ll
      pvlm = fgl.v[fgl.IDs[l]]
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
    maxl = ASCIIString("")
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

function maxIvWithoutID(ivs::Dict{ASCIIString, Float64}, l::AbstractString)
  max = 0
  maxl = ASCIIString("")
  for i in ivs
    if max < i[2] && i[1] != l;  max = i[2]; maxl = i[1]; end
  end
  return maxl
end

# binary tests to distinguish how to automatically add a landmark to the existing factor graph
function doAutoEvalTests(fgl::FactorGraph, ivs::Dict{ASCIIString, Float64}, maxl::ASCIIString, lmid::Int, lmindx::Int)
  maxAnyval = maxl != ASCIIString("") ? ivs[maxl] : 0.0
  # maxid = fgl.IDs[maxl]
  lmidSugg = lmid != -1 # a landmark ID has been suggested
  maxAnyExists = maxAnyval > 0.03 # there is notable intersection with previous landm
  lmIDExists = haskey(fgl.v, lmid) # suggested lmid already in fgl
  newlmindx = lmindx
  if lmIDExists
    lmSuggLbl = ASCIIString(fgl.v[lmid].label)
  else
    newlmindx = lmindx + 1
    lmSuggLbl = ASCIIString(string('l',newlmindx))
  end
  maxl2 = lmIDExists ? maxIvWithoutID(ivs, lmSuggLbl) : ASCIIString("")
  maxl2Exists = lmIDExists ? (maxl2 != "" ? ivs[maxl2] > 0.03 : false) : false # there is notable intersection with previous landm
  intgLmIDExists = lmIDExists ? ivs[lmSuggLbl] > 0.03 : false

  return lmidSugg, maxAnyExists, maxl2Exists, maxl2, lmIDExists, intgLmIDExists, lmSuggLbl, newlmindx
end

# TODO, busy here
function evalAutoCases!(fgl::FactorGraph, lmid::Int, ivs::Dict{ASCIIString, Float64}, maxl::ASCIIString,
                        pose::ASCIIString, lmPts::Array{Float64,2}, br::Array{Float64,1}, cov::Array{Float64,2}, lmindx::Int;
                        N::Int=100)
  lmidSugg, maxAnyExists, maxl2Exists, maxl2, lmIDExists, intgLmIDExists, lmSuggLbl, newlmindx = doAutoEvalTests(fgl,ivs,maxl,lmid, lmindx)

  println("evalAutoCases -- found=$(lmidSugg), $(maxAnyExists), $(maxl2Exists), $(lmIDExists), $(intgLmIDExists)")

  vlm = Union{}; fbr = Union{};
  if (!lmidSugg && !maxAnyExists)
    #new landmark and UniBR constraint
    v,L,lm = getLastLandm2D(fgl)
    vlm = newLandm!(fgl, lm, lmPts, cov, N=N)
    fbr = addBRFG!(fgl, pose, lm, br, cov)
  elseif !lmidSugg && maxAnyExists
    # add UniBR to best match maxl
    vlm = fgl.v[fgl.IDs[maxl]]
    fbr = addBRFG!(fgl, pose, maxl, br, cov)
  elseif lmidSugg && !maxl2Exists && !lmIDExists
    #add new landmark and add UniBR to suggested lmid
    vlm = newLandm!(fgl, lmSuggLbl, lmPts, cov, N=N)
    fbr = addBRFG!(fgl, pose, lmSuggLbl, br, cov)
  elseif lmidSugg && !maxl2Exists && lmIDExists && intgLmIDExists
    # doesn't self intesect with existing lmid, add UniBR to lmid
    vlm = fgl.v[lmid]
    fbr = addBRFG!(fgl, pose, lmSuggLbl, br, cov)
  elseif lmidSugg && maxl2Exists && !lmIDExists
    # add new landmark and add MMBR to both maxl and lmid
    vlm = newLandm!(fgl, lmSuggLbl, lmPts, cov, N=N)
    addMMBRFG!(fgl, pose, [maxl2;lmSuggLbl], br, cov)
  elseif lmidSugg && maxl2Exists && lmIDExists && intgLmIDExists
    # obvious case, add MMBR to both maxl and lmid. Double intersect might be the same thing
    println("evalAutoCases! -- obvious case is happening")
    addMMBRFG!(fgl, pose, [maxl2;lmSuggLbl], br, cov)
    vlm = fgl.v[fgl.IDs[lmSuggLbl]]
  elseif lmidSugg && maxl2Exists && lmIDExists && !intgLmIDExists
    # odd case, does not intersect with suggestion, but does with some previous landm
    # add MMBR
    warn("evalAutoCases! -- no self intersect with suggested $(lmSuggLbl) detected")
    addMMBRFG!(fgl, pose, [maxl;lmSuggLbl], br, cov)
    vlm = fgl.v[fgl.IDs[lmSuggLbl]]
  elseif lmidSugg && !maxl2Exists && lmIDExists && !intgLmIDExists
  #   # landm exists but no intersection with existing or suggested lmid
  #   # may suggest some error
    warn("evalAutoCases! -- no intersect with suggested $(lmSuggLbl) or map detected, adding  new landmark MM constraint incase")
    v,L,lm = getLastLandm2D(fgl)
    vlm = newLandm!(fgl, lm, lmPts, cov, N=N)
    addMMBRFG!(fgl, pose, [lm; lmSuggLbl], br, cov)
  else
    error("evalAutoCases! -- unknown case encountered, can reduce to this error to a warning and ignore user request")
  end

  return vlm, fbr, newlmindx
end

function addAutoLandmBR!(fgl::FactorGraph, pose::ASCIIString, lmid::Int, br::Array{Float64,1}, cov::Array{Float64,2}, lmindx::Int;
                      N::Int=100)
    vps = fgl.v[fgl.IDs[pose]]
    lmPts = projNewLandmPoints(vps, br, cov)
    lmkde = kde!(lmPts)
    currage = parse(Int, pose[2:end])
    ivs, maxl = calcIntersectVols(fgl, lmkde, currage=currage,maxdeltaage=6)

    # There are 8 cases of interest
    vlm, fbr, newlmindx = evalAutoCases!(fgl, lmid, ivs, maxl,pose,lmPts, br,cov,lmindx,N=N)

    # OLD TESTS
    # maxval = maxl != ASCIIString("") ? ivs[maxl] : 0.0
    # println("addAutoLandm! -- max intg val $(maxval)")
    # lm = Union{}; vlm = Union{};
    # if maxval > 0.03 # TODO : needs to be properly normalized.
    #   vlm = fgl.v[fgl.IDs[maxl]]
    #   lm = maxl
    #   println("prev lm age=$(vlm.attributes["maxage"])")
    #   # fbr = addBRFG!(fgl, pose, maxl, br, cov)
    #   # return vlm, fbr
    # else
    #   v,L,lm = getLastLandm2D(fgl)
    #   vlm = newLandm!(fgl, lm, lmPts, cov, N=N)
    # end
    # fbr = addBRFG!(fgl, pose, lm, br, cov)
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

function initFactorGraph!(fg::FactorGraph; P0=diagm([0.03;0.03;0.001]),
                      init=[0.0;0.0;0.0], N::Int=100, lbl="x1")
    init = (init')'
    v1 = addNode!(fg, lbl, init, P0, N=N)
    addFactor!(fg, PriorPose2([v1], init, P0,  [1.0]) )
    return lbl
end
