# utility functions that provide Gaussian odometry accumulation

import IncrementalInference: getFactorMean

export getFactorMean
export accumulateDiscreteLocalFrame!, duplicateToStandardFactorVariable, extractDeltaOdo, resetFactor!
export odomKDE


getFactorMean(fct::PriorPose2) = getFactorMean(fct.Z)
getFactorMean(fct::Pose2Pose2) = getFactorMean(fct.z)
getFactorMean(fct::MutablePose2Pose2Gaussian) = getFactorMean(fct.Zij)


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
                                       dt::Float64=1.0;
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

Helper function to duplicate values from a special factor variable into standard factor and variable.  Returns the name of the new factor.

Notes:
- Developed for accumulating odometry in a `MutablePosePose` and then cloning out a standard PosePose and new variable.
- Does not change the original MutablePosePose source factor or variable in any way.
- Assumes timestampe from mpp object.
"""
function duplicateToStandardFactorVariable(::Type{Pose2Pose2},
                                           mpp::MutablePose2Pose2Gaussian,
                                           dfg::AbstractDFG,
                                           prevsym::Symbol,
                                           newsym::Symbol;
                                           solvable::Int=1,
                                           autoinit::Union{Nothing, Bool}=nothing,
                                           graphinit::Bool=true,
                                           cov::Union{Nothing, Matrix{Float64}}=nothing  )::Symbol
  #
  graphinit = autoinit != nothing ? (@warn "autoinit deprecated, use graphinit instead."; autoinit) : graphinit

  # extract factor values and create PosePose object
  posepose = Pose2Pose2(MvNormal(mpp.Zij.μ, cov==nothing ? mpp.Zij.Σ.mat : cov))

  # modify the factor graph
  addVariable!(dfg, newsym, Pose2, solvable=solvable, timestamp=mpp.timestamp)
  addFactor!(dfg, [prevsym; newsym], posepose, solvable=solvable, graphinit=graphinit, timestamp=mpp.timestamp)
  return ls(dfg, newsym)[1]
end

"""
    $SIGNATURES

Reset the transform value stored in a `::MutablePose2Pose2Gaussian` to zero.
"""
function resetFactor!(mpp::MutablePose2Pose2Gaussian)::Nothing
  mpp.Zij = MvNormal(zeros(3), 1e-6*Matrix{Float64}(LinearAlgebra.I, 3,3) )
  nothing
end


"""
    $SIGNATURES

Extract deltas from existing dead reckoning data so that odometry calculations can be repeated later.

Notes
- Useful for reverse engineering data or simulation tools.
"""
function extractDeltaOdo(XX, YY, TH)
  dt = 1.0
  DX = zeros(3,length(XX))
  nXYT__ = zeros(3,size(DX,2))
  nXYT__[:,1] = [XX[1];YY[1];TH[1]]
  for i in 2:length(XX)
    wTbk = SE2([XX[i-1];YY[i-1];TH[i-1]])
    wTbk1 = SE2([XX[i];YY[i];TH[i]])
    bkTbk1 = wTbk\wTbk1
    DX[:,i] = se2vee(bkTbk1)

    # test
    nXYT__[:,i] .= se2vee(SE2(nXYT__[:,i-1])*SE2(DX[:,i]))
    # nXYT__[:,i] .= se2vee(SE2(nXYT__[:,i-1])*bkTbk1)
  end

  return DX
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
