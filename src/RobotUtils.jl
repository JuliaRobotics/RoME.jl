
include("ISAMRemoteSolve.jl")

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

function getLastPose2D(fg::FactorGraph)
  max = 0
  maxid = 0
  for v in fg.v
      if v[2].attributes["label"][1] == 'x'
        val = parse(Int,v[2].attributes["label"][2:end])
        if max < val
          max = val
          maxid = v[1]
        end
      end
  end
  v = fg.v[maxid]
  X = getVal(v)
  return v, X, string('x',max+1)
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

# cov should not be required here
function newLandm!(fg::FactorGraph, lm::ASCIIString, wPos::Array{Float64,2}, sig::Array{Float64,2};
                  N::Int=100)
    v=addNode!(fg, lm, wPos, sig, N=N)
    return v
end

function addBRFG!(fg::FactorGraph, pose::ASCIIString,
                  lm::ASCIIString, br::Array{Float64,1}, cov::Array{Float64,2})
    vps = fg.v[fg.IDs[pose]]
    vlm = fg.v[fg.IDs[lm]]
    f = addFactor!(fg, Pose2DPoint2DBearingRange([vps;vlm],(br')',  cov,  [1.0]) )
    u, P = pol2cart(br[[2;1]], diag(cov))
    infor = inv(P^2)
    addLandmMeasRemote(vps.index,vlm.index,u,infor)
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

function projNewLandm!(fg::FactorGraph, pose::ASCIIString, lm::ASCIIString, br::Array{Float64,1}, cov::Array{Float64,2};
                        addfactor=true, N::Int=100)
    vps = fg.v[fg.IDs[pose]]
    Xps = getVal(vps)
    lmPts = zeros(2,size(Xps,2))
    for i in 1:size(Xps,2)
        ent = [cov[1,1]*randn(); cov[2,2]*randn()]
        init = vec(Xps[1:2,i])+randn(2)
        lmPts[:,i] = solveLandm(br + ent, vec(Xps[:,i]), init)
    end
    vlm = newLandm!(fg, lm, lmPts, cov, N=N) # cov should not be required here
    if addfactor
      fbr = addBRFG!(fg, pose, lm, br, cov)
      return vlm, fbr
    else
      return vlm
    end
end

function malahanobisBR(measA, preA, cov::Array{Float64,2})
    # measure landmark with noise
    res = measA - preA
    lambda = cov \ eye(2)
    #Malahanobis distance
    mala2 = res' * lambda * res
    mala = sqrt(mala2)
    return mala
end

function initFactorGraph!(fg::FactorGraph; P0=diagm([0.01;0.01;0.001]),
                      init=[0.0;0.0;0.0], N::Int=100, lbl="x1")
    init = (init')'
    v1 = addNode!(fg, lbl, init, P0, N=N)
    addFactor!(fg, PriorPose2([v1], init, P0,  [1.0]) )
    return lbl
end
