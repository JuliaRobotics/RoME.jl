# utility functions that provide Gaussian odometry accumulation

export accumulateDiscreteLocalFrame!
export rebaseFactorVariable!, duplicateToStandardFactorVariable
export solveBinaryFactorParameteric, getFactorMean, accumulateFactorMeans


function cont2disc(F::Matrix{Float64},
                   G::Matrix{Float64},
                   Qc::Matrix{Float64},
                   dt::Float64,
                   Phik::Matrix{Float64}=Matrix{Float64}(LinearAlgebra.I, 0,0) )
    #
    fr,fc = size(F)
    gr,gc = size(G)

    # two bits of new memory allocated
    M1 = zeros(fc+gc,fc+gc)
    M2 = zeros(fr+fc,fr+fc)

    M1[1:fr,1:fc] = F
    M1[1:gr,(fc+1):end] = G #1:gr,(fc+1):(fc+gc)

    # must convert to propagateLinSystem call, use trapezoidal
    Md1 = exp(M1*dt) # heavy lifting here
    Phi = size(Phik,1) == 0 ? Md1[1:fr,1:fc] : Phik
    Gamma = Md1[1:fr,(fc+1):end]

    #M2 = [[-F';(G*Qc*G')']';[zeros(9,9);F]'] # easy concat
    GQG = (G*Qc*G')
    gqgr, gqgc = size(GQG)
    M2[1:fc,1:fr] = -F
    M2[1:fr,(fc+1):end] = GQG
    M2[(fr+1):end,(fc+1):end] = F'

    Md2 = exp(M2*dt) # heavy lifting here
    Qd = Phi * Md2[1:fr,(fc+1):end] #Qd = Phi*(Md2[1:fr,(fc+1):end])

    # Qd = GQG*dt;

    return Phi, Gamma, Qd
end


"""
    $SIGNATURES

Advance an odometry factor as though integrating an ODE -- i.e. ``X_2 = X_1 ⊕ ΔX``. Accepts continuous domain process noise density `Qc` which is internally integrated to discrete process noise Qd.  ``DX`` is assumed to already be incrementally integrated before this function.  See related `accumulateContinuousLocalFrame!` for fully continuous system propagation.

Notes
- This update stays in the same reference frame but updates the local vector as though accumulating measurement values over time.
- Kalman filter would have used for noise propagation: ``Pk1 = F*Pk*F' + Qdk``
- From Chirikjian, Vol.II, 2012, p.35: Jacobian SE(2), Jr = [cθ sθ 0; -sθ cθ 0; 0 0 1] -- i.e. dSE2/dX' = SE2([0;0;-θ])
- `DX = dX/dt*Dt`
- assumed process noise for `{}^b Qc = {}^b [x;y;yaw] = [fwd; sideways; rotation.rate]`

Dev Notes
- TODO many operations here can be done in-place.

Related

accumulateContinuousLocalFrame!, accumulateDiscreteReferenceFrame!
"""
function accumulateDiscreteLocalFrame!(mpp::MutablePose2Pose2Gaussian,
                                       DX::Vector{Float64},
                                       Qc::Matrix{Float64},
                                       dt::Float64;
                                       Fk = SE2([0;0;-DX[3]]),
                                       Gk = Matrix{Float64}(LinearAlgebra.I, 3,3),
                                       Phik = SE2(DX) )
  #
  kXk1 = SE2(mpp.Zij.μ)*Phik
  phi, gamma, Qd = cont2disc(Fk, Gk, Qc, dt, Phik)
  Covk1 = Phik*(mpp.Zij.Σ.mat)*(Phik') + Qd
  check = norm(Covk1 - Covk1')
  1e-4 < check ? @warn("Covk1 is getting dangerously non-Hermitian, still forcing symmetric covariance matrix.") : nothing
  @assert check < 1.0
  Covk1 .+= Covk1'
  Covk1 ./= 2
  mpp.Zij = MvNormal(se2vee(kXk1), Covk1)
  nothing
end

"""
    $SIGNATURES

Helper function to duplicate values from a special factor variable into standard factor and variable.

Notes:
- Developed for accumulating odometry in a `MutablePosePose` and then cloning out a standard PosePose and new variable.
- Does not change the original MutablePosePose source factor or variable in any way.
"""
function duplicateToStandardFactorVariable(::Pose2Pose2,
                                           mpp::MutablePose2Pose2Gaussian,
                                           dfg::AbstractDFG,
                                           prevsym::Symbol,
                                           newsym::Symbol )::Nothing
  #
  # extract factor values and create PosePose object
  posepose = Pose2Pose2(deepcopy(mpp.Zij))

  # modify the factor graph
  addVariable!(dfg, newsym, Pose2)
  addFactor!(dfg, [prevsym; newsym], posepose)
  nothing
end



"""
    $SIGNATURES

Helper function to modify factor connectivity to variables.

Notes
- Developed for updating a dead reckoning odometry factor.
- Arguments are order sensitive.
"""
function rebaseFactorVariable!(dfg::AbstractDFG,
                               fctsym::Symbol,
                               newvars::Vector{Symbol};
                               rmDisconnected::Bool=true,
                               autoinit::Bool=false  )::Nothing
  #
  # check that all new variables are available
  @assert sum(map(x->hasVariable(dfg, x), newvars)) == length(newvars)

  # get existing factor details
  fct = getFactor(dfg, fctsym)
  fcttype = getFactorType(fct)
  mh = getMultihypoDistribution(fct)

  # get old vars
  oldvars = getVariableOrder(fct)

  # delete old factor from graph
  deleteFactor!(dfg, fctsym)

  # add the factor back into graph against new variables
  addFactor!(dfg, newvars, fcttype, autoinit=autoinit, multihypo=mh)

  # clean up disconnected variables if requested
  if rmDisconnected
    for ov in oldvars
      # find variables that are not connected to anything
      if length(ls(dfg, ov)) == 0
        deleteVariable!(dfg, ov)
      end
    end
  end

  return nothing
end


## Previous methods

function odomKDE(p1,dx,cov)
  @warn "odomKDE is beig deprecated in its current form, consider using approxConv or predictVariableByFactor instead."
  X = getPoints(p1)
  sig = diag(cov)
  RES = zeros(size(X))
  # increases the number of particles based on the number of modes in the measurement Z
  for i in 1:size(X,2)
      ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
      RES[:,i] = addPose2Pose2(X[:,i], dx + ent)
  end
  return manikde!(RES, Pose2)
end



function getFactorMean(fct::FunctorInferenceType)
  fctt = typeof(getFactorType(fct))
  error("no getFactorMean defined for $(fctt.name), has fields $(fieldnames(fctt))")
end

getFactorMean(fct::MvNormal) = fct.μ
getFactorMean(fct::BallTreeDensity) = getKDEMean(fct)
getFactorMean(fct::AliasingScalarSampler) = Statistics.mean(rand(fct,1000))

getFactorMean(fct::PriorPose2) = getFactorMean(fct.Z)
getFactorMean(fct::Pose2Pose2) = getFactorMean(fct.z)
getFactorMean(fct::MutablePose2Pose2Gaussian) = getFactorMean(fct.Zij)


getFactorMean(fct::DFGFactor) = getFactorMean(getFactorType(fct))

getFactorMean(dfg::AbstractDFG, fctsym::Symbol) = getFactorMean(getFactor(dfg, fctsym))

"""
    $SIGNATURES

Helper function to propagate a parametric estimate along a factor chain.

Notes
- Not used during mmisam inference.
- Expected uses are for user analysis of factors and estimates.
- real-time dead reckoning chain prediction.

Related:

approxConv, accumulateFactorMeans, MutablePose2Pose2Gaussian
"""
function solveBinaryFactorParameteric(dfg::AbstractDFG,
                                      fct::DFGFactor,
                                      currval::Vector{Float64},
                                      srcsym::Symbol,
                                      trgsym::Symbol  )::Vector{Float64}
  #
  outdims = getVariableDim(getVariable(dfg, trgsym))
  meas = getFactorType(fct)
  mea = getFactorMean(fct)
  measT = (reshape(mea,:,1),)

  # calculate the projection
  varmask = (1:2)[getVariableOrder(fct) .== trgsym][1]
  pts = approxConvBinary( reshape(currval,:,1), meas, outdims, measT, varidx=varmask )

  # return the result
  @assert length(pts) == outdims
  return pts[:]
end

"""
    $SIGNATURES

Accumulate chains of binary factors---potentially starting from a prior---as a parameteric mean value only.

Notes
- Not used during mmisam inference.
- Expected uses are for user analysis of factors and estimates.
- real-time dead reckoning chain prediction.

Related:

approxConv, solveBinaryFactorParameteric, MutablePose2Pose2Gaussian
"""
function accumulateFactorMeans(dfg::AbstractDFG, fctsyms::Vector{Symbol})

  ## get the starting estimate
  val = zeros(0)
  nextidx = 1
  onePrior = false
  currsym = :null
  if isPrior(dfg, fctsyms[nextidx])
    # if first factor is prior
    @assert !onePrior
    onePrior = true
    val = getFactorMean(dfg, fctsyms[nextidx])
    currsym = ls(dfg, fctsyms[nextidx])
    nextidx += 1
  else
    # get first value from current variable estimate
    vars = getVariableOrder(dfg, fctsyms[nextidx])
    val = calcVariablePPE(dfg, vars[nextidx]).suggested
    currsym = setdiff(vars, intersect(vars, ls(dfg, fctsyms[nextidx+1])))[1]
  end

  # Propagate the parametric value along the factor chain
  for fct in map(x->getFactor(dfg, x), fctsyms[nextidx:end])
    # first find direction of solve
    vars = getVariableOrder(fct)
    srcsym = currsym
    trgsym = setdiff(vars, [srcsym])[1]
    varmask = (1:2)[getVariableOrder(fct) .== trgsym][1]
    val = solveBinaryFactorParameteric(dfg,fct,val,:x0,:x1)
  end

  return val
end
