# Bearing and Range constraints for 2D


# better to use bearingrange with [uniform bearing], numerical solving issue on 1D
type Pose2DPoint2DRange <: IncrementalInference.FunctorPairwise
    Zij::Vector{Float64} # bearing and range hypotheses as columns
    Cov::Float64
    W::Vector{Float64}
    Pose2DPoint2DRange() = new()
    Pose2DPoint2DRange(x...) = new(x[1],x[2],x[3])
end
function (pp2r::Pose2DPoint2DRange)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}, Array{Float64,1}}, # from getSample
      xi::Array{Float64,2},
      lm::Array{Float64,2}  )
  #
  # TODO -- still need to add multi-hypotheses support here
  # this is the noisy range
  z = pp2r.Zij[1]+meas[1][1,idx]
  # theta = meas[2]
  # @show size(lm), size(xi), size(meas), size(meas[1]), size(meas[2])
  XX = lm[1,idx] - (z*cos(meas[2][idx]) + xi[1,idx])
  YY = lm[2,idx] - (z*sin(meas[2][idx]) + xi[2,idx])
  res[1] = XX^2 + YY^2
  nothing
end
function getSample(pp2::Pose2DPoint2DRange, N::Int=1)
  return (pp2.Cov*randn(1,N),  2*pi*rand(N))
end


#-------------------------------------------------------------------------------
# bearing and range available

type Pose2DPoint2DBearingRange{B <: Distributions.Distribution, R <: Distributions.Distribution} <: IncrementalInference.FunctorPairwise
    # Zij::Array{Float64,2} # bearing and range hypotheses as columns
    # Cov::Array{Float64,2}
    # W::Array{Float64,1}
    bearing::B
    range::R
    Pose2DPoint2DBearingRange() = new()
    Pose2DPoint2DBearingRange{B,R}(x1::B,x2::R) = new(x1,x2)
end
function getSample(pp2br::Pose2DPoint2DBearingRange, N::Int=1)
  b = rand(pp2br.bearing, N)
  r = rand(pp2br.range, N)
  return ([b';r'], )
end
# define the conditional probability constraint
function (pp2br::Pose2DPoint2DBearingRange)(res::Array{Float64},
        idx::Int,
        meas::Tuple{Array{Float64,2}},
        xi::Array{Float64,2},
        lm::Array{Float64,2} )
  #
  # @show size(lm), size(xi), size(meas), size(meas[1]), size(meas[2])
  res[1] = lm[1,idx] - (meas[1][2,idx]*cos(meas[1][1,idx]) + xi[1,idx])
  res[2] = lm[2,idx] - (meas[1][2,idx]*sin(meas[1][1,idx]) + xi[2,idx])
  nothing
end




# Support for database based solving

passTypeThrough(d::FunctionNodeData{Pose2DPoint2DRange}) = d

type PackedPose2DPoint2DBearingRange <: IncrementalInference.PackedInferenceType
    bmu::Float64 # 0rotations, 1translation in each column
    bsig::Float64
    rmu::Float64
    rsig::Float64
    PackedPose2DPoint2DBearingRange() = new()
    PackedPose2DPoint2DBearingRange(x...) = new(x[1], x[2], x[3], x[4])
end
function convert(::Type{Pose2DPoint2DBearingRange{Normal, Normal}}, d::PackedPose2DPoint2DBearingRange)
  return Pose2DPoint2DBearingRange{Distributions.Normal, Distributions.Normal}(Distributions.Normal(d.bmu,d.bsig), Distributions.Normal(d.rmu, d.rsig))
end
function convert(::Type{PackedPose2DPoint2DBearingRange}, d::Pose2DPoint2DBearingRange{Normal, Normal})
  return PackedPose2DPoint2DBearingRange(d.bearing.μ, d.bearing.σ,
                                         d.range.μ,   d.range.σ )
end



# to be deprecated
function pack3(xL1, xL2, p1, p2, p3, xF3)
    warn("RoME.BearingRange2D:pack3 to be deprecated")
    X = zeros(3)
    X[p1] = xL1
    X[p2] = xL2
    X[p3] = xF3
    return X
end

function bearrang!(residual::Array{Float64,1}, Z::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
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
    p = collect(1:3);
    shuffle!(p);
    p1 = p.==1; p2 = p.==2; p3 = p.==3
    #@show init, par
    r = nlsolve(    (x, res) -> bearrang!(res, Zbr,  pack3(x[1], x[2], p1, p2, p3, init[p3]), par),
                    [init[p1];init[p2]] )
    return pack3(r.zero[1], r.zero[2], p1, p2, p3, init[p3]);
end

function solveSetSeps(fnc::Function, Zbr::Array{Float64,1}, CovZ::Array{Float64,2},
                      pars::Array{Float64,2}, inits::Array{Float64,2})
    out = zeros(size(inits))
    for i in 1:size(pars,2)
        ent = 1.0*[CovZ[1,1]*randn(); CovZ[2,2]*randn()]
        out[:,i] = fnc((Zbr+ent), vec(pars[:,i]), vec(inits[:,i]) )
    end
    return out
end

# Xid is the one you want to get back
function evalPotential(brpho::Pose2DPoint2DBearingRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
    # TODO -- add null hypothesis here, might even be done one layer higher in call stack
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

function evalPotentialNew(brpho::Pose2DPoint2DBearingRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
    # TODO -- add null hypothesis here, might even be done one layer higher in call stack
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
function evalPotential(rho::Pose2DPoint2DRange, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
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




# no longer needed
# function convert(::Type{PackedFunctionNodeData{PackedPose2DPoint2DBearingRange}}, d::FunctionNodeData{Pose2DPoint2DBearingRange})
#   return PackedFunctionNodeData{PackedPose2DPoint2DBearingRange}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           string(d.frommodule), convert(PackedPose2DPoint2DBearingRange, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Pose2DPoint2DBearingRange}}, d::PackedFunctionNodeData{PackedPose2DPoint2DBearingRange})
#   return FunctionNodeData{Pose2DPoint2DBearingRange}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           Symbol(d.frommodule), convert(Pose2DPoint2DBearingRange, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{Pose2DPoint2DBearingRange})
#   return convert(PackedFunctionNodeData{PackedPose2DPoint2DBearingRange}, d)
# end
# function FNDdecode(d::PackedFunctionNodeData{PackedPose2DPoint2DBearingRange})
#   return convert(FunctionNodeData{Pose2DPoint2DBearingRange}, d)
# end


# ------------------------------------------------------
