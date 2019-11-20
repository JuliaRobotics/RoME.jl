# utility functions that provide Gaussian odometry accumulation

export accumulateDiscreteLocalFrame!
export rebaseFactorVariable!, duplicateToStandardFactorVariable



"""
    $SIGNATURES

Advance an odometry factor as though integrating an ODE -- i.e. ``X_2 = X_1 ⊕ ΔX``.

Notes
- This update stays in the same reference frame but updates the local vector as though accumulating measurement values over time.
- Kalman filter would have used for noise propagation: ``Pk1 = F*Pk*F' + Qdk``

Related

accumulateDiscreteReferenceFrame!
"""
function accumulateDiscreteLocalFrame!(mpp::MutablePose2Pose2Gaussian,
                                       DX::Vector{Float64},
                                       Qd::Matrix{Float64})
  #
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
                               newvars::Vector{Symbols};
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
