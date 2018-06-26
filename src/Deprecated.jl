
# old code below

# This has been moved to transform utils
# TODO Switch to using SE(2) oplus
# DX = [transx, transy, theta]
function addPose2Pose2!(retval::Array{Float64,1}, x::Array{Float64,1}, dx::Array{Float64,1})
  X = SE2(x)
  DX = SE2(dx)
  se2vee!(retval, X*DX)
  nothing
end
function addPose2Pose2(x::Array{Float64,1}, dx::Array{Float64,1})
    retval = zeros(3)
    addPose2Pose2!(retval, x, dx)
    return retval
end


function evalPotential(obs::PriorPose2, Xi::Array{Graphs.ExVertex,1}; N::Int=200)
    cov = diag(obs.Cov)
    ret = zeros(3,N)
    warn("should not be running")
    for j in 1:N
      for i in 1:size(obs.Zi,1)
        ret[i,j] += obs.Zi[i,1] + (cov[i]*randn())
      end
    end
    return ret
end



function evalPotential(odom::Pose2Pose2, Xi::Array{Graphs.ExVertex,1}, Xid::Int; N::Int=100)
    rz,cz = size(odom.Zij)
    Xval = Array{Float64,2}()
    XvalNull = Array{Float64,2}()
    warn("should not be running")
    # implicit equation portion -- bi-directional pairwise function
    if Xid == Xi[1].index #odom.
        #Z = (odom.Zij\eye(rz)) # this will be used for group operations
        Z = se2vee(SE2(vec(odom.Zij)) \ eye(3))
        Xval = getVal(Xi[2])
        XvalNull = getVal(Xi[1])
    elseif Xid == Xi[2].index
        Z = odom.Zij
        Xval = getVal(Xi[1])
        XvalNull = getVal(Xi[2])
    else
        error("Bad evalPairwise Pose2Pose2")
    end

    r,c = size(Xval)
    RES = zeros(r,c*cz)

    # TODO -- this should be the covariate error from Distributions, only using diagonals here (approxmition for speed in first implementation)
    # dd = size(Z,1) == r
    ENT = randn(r,c)
    HYP = rand(Categorical(odom.W),c) # TODO consolidate
    HYP -= length(odom.W)>1 ? 1 : 0
    for d in 1:r
       @fastmath @inbounds ENT[d,:] = ENT[d,:].*odom.Cov[d,d]
    end
    # Repeat ENT values for new modes from meas
    for j in 1:cz
      for i in 1:c
        if HYP[i]==1 # TODO consolidate hypotheses on Categorical
          z = Z[1:r,j].+ENT[1:r,i]
          RES[1:r,i*j] = addPose2Pose2(Xval[1:r,i], z )
        else
          RES[1:r,i*j] = XvalNull[1:r,i]
        end
      end
    end

    return RES
end



# to be deprecated
function pack3(xL1, xL2, p1, p2, p3, xF3)
    error("RoME.BearingRange2D:pack3 to be deprecated")
    X = zeros(3)
    X[p1] = xL1
    X[p2] = xL2
    X[p3] = xF3
    return X
end

function bearrang!(residual::Array{Float64,1}, Z::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
  warn("bearrang! is deprecated")
  wTb = SE2(X)
  bTl = wTb\[L[1:2];1.0]
  b = atan2(bTl[2],bTl[1])
  residual[1] = Z[2]-norm(bTl[1:2])
  residual[2] = Z[1]-b
  nothing
end


# should be collapsed to use only numericRootGenericRandomizedFnc
function solveLandm(Zbr::Array{Float64,1}, par::Array{Float64,1}, init::Array{Float64,1})
    return numericRoot(bearrang!, Zbr, par, init)
    # return (nlsolve(   (l, res) -> bearrang!(res, Zbr, par, l), init )).zero
end

# old numeric residual function for pose 2 to pose 2 constraint function.
function solvePose2(Zbr::Array{Float64,1}, par::Array{Float64,1}, init::Array{Float64,1})
    # TODO -- rework to ominus oplus and residual type method
    error("solvePose2 is deprecated")
    p = collect(1:3);
    shuffle!(p);
    p1 = p.==1; p2 = p.==2; p3 = p.==3
    #@show init, par
    r = nlsolve(    (res, x) -> bearrang!(res, Zbr,  pack3(x[1], x[2], p1, p2, p3, init[p3]), par),
                    [init[p1];init[p2]] )
    return pack3(r.zero[1], r.zero[2], p1, p2, p3, init[p3]);
end

function solveSetSeps(fnc::Function, Zbr::Array{Float64,1}, CovZ::Array{Float64,2},
                      pars::Array{Float64,2}, inits::Array{Float64,2})
    error("solveSetSeps is deprecated")
    out = zeros(size(inits))
    for i in 1:size(pars,2)
        ent = 1.0*[CovZ[1,1]*randn(); CovZ[2,2]*randn()]
        out[:,i] = fnc((Zbr+ent), vec(pars[:,i]), vec(inits[:,i]) )
    end
    return out
end

# Xid is the one you want to get back
function evalPotential(brpho::Pose2DPoint2DBearingRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int; N::Int=100)
    # TODO -- add null hypothesis here, might even be done one layer higher in call stack
    error("evalPotential(brpho::Pose2DPoint2DBearingRange,...) should not be here anymore")
    val = Array{Float64,2}()
    ini = Array{Graphs.ExVertex,1}()
    par = Array{Graphs.ExVertex,1}()
    oth = Array{Graphs.ExVertex,1}()
    ff::Function = +
    nullhyp = 0.0
    # mmodes = 1 < length(brpho.W)
    # implicit equation portion -- multi-dependent function
    if Xid == Xi[1].index # brpho. ## find the pose

        ff = solvePose2
        # ini = brpho.Xi[1]
        par = Xi[2:end]
        for j in 1:(length(Xi)-1)
            push!(ini, Xi[1])
        end
        #println("Xid == brpho.Xi[1].index inits=", size(inits), ", par=", size(pars))
    elseif Xid == Xi[2].index # find landmark
        if length(Xi) > 2
            nullhyp = 0.5 # should be variable weight
            oth = Xi[3]
        end
        ff = solveLandm
        ini = Xi[2]
        par = Xi[1]
    elseif length(Xi) > 2
        nullhyp = 0.5 # should be variable weight
        if Xid == Xi[3].index # find second mode landmark
            ff = solveLandm
            ini = Xi[3]
            oth = Xi[2]
            par = Xi[1]
        end
    end
    if ff == +
        error("Bad evalPotential Pose2Point2DBearingRange")
    end

    # Gamma = Categorical(brpho.W)

    inits = getVal(ini)
    pars = getVal(par)
    others =  length(Xi) > 2 && Xid != Xi[1].index ? getVal(oth) : Union{}
    # add null hypothesis case
    len = length(Xi) > 2 && Xid != Xi[1].index ? size(others,2) : 0
    # gamma = mmodes ? rand(Gamma) : 1
    numnh = floor(Int, 2*nullhyp*len) # this doubles the value count for null cases
    nhvals = zeros(size(inits,1),numnh)
    for i in 1:numnh
        idx = floor(Int,len*rand()+1)
        # nhvals[:,i] = inits[:,idx] # WRONG! this should be the other landmark, not current landmark
        nhvals[:,i] = others[:,idx]
    end

    val = solveSetSeps(ff, vec(brpho.Zij[:,1]), brpho.Cov, pars, inits)

    return [val';nhvals']'
end

function evalPotentialNew(brpho::Pose2DPoint2DBearingRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int; N::Int=100)
    # TODO -- add null hypothesis here, might even be done one layer higher in call stack
    error("deprecated")
    val = Array{Float64,2}()
    ini = Array{Graphs.ExVertex,1}()
    par = Array{Graphs.ExVertex,1}()
    ff::Function = +
    mmodes = 1 < length(brpho.W)
    x = getVal(brpho.Xi[1])
    l1 = getVal(brpho.Xi[2])
    l2 = mmodes ? getVal(brpho.Xi[3]) : nothing
    nPts = size(x,2)

    pars = Array{Float64,2}(size(x))

    # discrete option, marginalized out before message is sent
    Gamma = mmodes ? rand(Categorical(brpho.W),nPts) : ones(Int,nPts)
    L1s = Gamma .== 1
    nl1s = sum(map(Int, L1s))
    nl2s = nPts-nl1s
    # implicit equation portion -- multi-dependent function

    if Xid == brpho.Xi[1].index # find the pose
        ff = solvePose2
        # par = brpho.Xi[2:end]
        push!(ini, brpho.Xi[1])
        for j in 1:nPts
            pars[:,j] = L1s[j] ? l1[:,j] : l2[:,j]
        end
    elseif Xid == brpho.Xi[2].index # find landmark
        ff = solveLandm
        pars = x
        push!(ini,brpho.Xi[2])
    elseif mmodes
        if Xid == brpho.Xi[3].index # find second mode landmark
            ff = solveLandm
            push!(ini,brpho.Xi[3])
            pars = x
        end
    end
    if ff == +
        error("Bad evalPotential Pose2Point2DBearingRange")
    end

    inits = getVal(ini)

    L1 = sample(getKDE(brpho.Xi[2]),numsampls) #?? # TODO you where here

    len = size(inits,2)
    numnh = floor(Int, 2*nullhyp*len) # this doubles the value count for null cases
    nhvals = zeros(size(inits,1),numnh)
    for i in 1:numnh
        idx = floor(Int,len*rand()+1)
        nhvals[:,i] = inits[:,idx]
    end

    val = solveSetSeps(ff, vec(brpho.Zij[:,1]), brpho.Cov, pars, inits)

    return [val';nhvals']'
end



# Solve for Xid, given values from vertices [Xi] and measurement rho
function evalPotential(rho::Pose2DPoint2DRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int; N::Int=100)
  error("deprecated")
  fromX, ret = nothing, nothing
  if Xi[1].index == Xid
    fromX = getVal( Xi[2] )
    ret = deepcopy(getVal( Xi[1] )) # carry pose yaw row over if required
    ret[3,:] = 2*pi*rand(size(fromX,2))-pi
  elseif Xi[2].index == Xid
    fromX = getVal( Xi[1] )
    ret = deepcopy(getVal( Xi[2] )) # carry pose yaw row over if required
  end
  r,c = size(fromX)
  theta = 2*pi*rand(c)
  noisy = rho.Cov*randn(c) + rho.Zij[1]

  for i in 1:c
    ret[1,i] = noisy[i]*cos(theta[i]) + fromX[1,i]
    ret[2,i] = noisy[i]*sin(theta[i]) + fromX[2,i]
  end

  return ret
end
