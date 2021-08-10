# canonical factor graph examples useful for development and learning.

export generateCanonicalFG_ZeroPose, generateCanonicalFG_Hexagonal, generateCanonicalFG_TwoPoseOdo, generateCanonicalFG_Circle
export warmUpSolverJIT

"""
    $SIGNATURES

Lower level utility function for simulated canonical graph generation, where a pose variable and factor is added using this function.

Notes
- Also adds the `refKey=:simulated` (default) PPE as in graph record of the simulated values. 
- Use callback `postpose_cb(g::AbstractDFG,lastpose::Symbol)` to call user operations after each pose step.

Related

[`generateCanonicalFG_ZeroPose`](@ref), [`generateCanonicalFG_Helix2D!`](@ref)
"""
function _addPoseCanonical!(fg::AbstractDFG, 
                            prevLabel::Symbol,
                            posecount::Int, # can be overriden with genLabel
                            factor::AbstractFactor;
                            poseRegex::Regex=r"x\d+",
                            genLabel = Symbol(match(r"[A-Za-z]+", poseRegex.pattern).match, posecount),
                            srcType::Type{<:InferenceVariable} = getVariableType(fg, prevLabel) |> typeof,
                            poseType::Type{<:InferenceVariable} = srcType, # control destination type TODO simplify
                            graphinit::Bool=false,
                            solvable::Integer=1,
                            variableTags::AbstractVector{Symbol}=Symbol[],
                            factorTags::AbstractVector{Symbol}=Symbol[],
                            refKey::Symbol=:simulated,
                            overridePPE=nothing,
                            postpose_cb::Function=(fg_,latestpose)->()  )
  #
  # calculate and add the reference value
  isAlready, simPPE, = IIF._checkVariableByReference(fg, prevLabel, poseRegex, poseType, factor, refKey=refKey, srcType=srcType, overridePPE=overridePPE)

  # dispatch on prior or binary factor
  _getlabels(fact::AbstractPrior) = [genLabel;]
  _getlabels(fact::AbstractRelative) = [prevLabel; genLabel]
  
  # add new pose variable
  v_n = addVariable!(fg, genLabel, poseType, solvable=solvable, tags=variableTags )
  addFactor!(fg, _getlabels(factor), factor, graphinit=graphinit, solvable=solvable, tags=factorTags )

  # store simulated PPE for future use
  # ppe = DFG.MeanMaxPPE(refKey, simPPE, simPPE, simPPE)
  setPPE!(v_n, refKey, typeof(simPPE), simPPE)

  # user callback in case something more needs to be passed down
  postpose_cb(fg, genLabel)

  # return the new variable
  return v_n
end




"""
    $SIGNATURES

Generate a canonical factor graph with a Pose2 `:x0` and MvNormal with covariance `P0`.

Notes
- Use e.g. `varType=Point2` to change from the default variable type `Pose2`.
- Use `priorArgs::Tuple` to override the default input arguments to `priorType`.
- Use callback `postpose_cb(g::AbstractDFG,lastpose::Symbol)` to call user operations after each pose step.
"""
function generateCanonicalFG_ZeroPose(; varType::Type{<:InferenceVariable}=Pose2,
                                        graphinit::Bool=true,
                                        dfg::AbstractDFG = LightDFG{SolverParams}(solverParams=SolverParams(graphinit=graphinit)),  
                                        useMsgLikelihoods::Bool=getSolverParams(dfg).useMsgLikelihoods,
                                        label::Symbol=:x0,
                                        priorType::Type{<:AbstractPrior}=DFG._getPriorType(varType),
                                        μ0::AbstractVector{<:Real}= zeros(getDimension(varType)),
                                        Σ0::AbstractMatrix{<:Real}= diagm(0.01*ones(getDimension(varType))),
                                        priorArgs::Tuple = (MvNormal(μ0, Σ0),),
                                        solvable::Integer=1,
                                        variableTags::AbstractVector{Symbol}=Symbol[],
                                        factorTags::AbstractVector{Symbol}=Symbol[],
                                        postpose_cb::Function=(fg_,latestpose)->()  )
  #

  # only add the first variable if none others exist
  if !exists(dfg, label)
    # generate a default prior
    prpo = priorType(priorArgs...)
    # add the variable and prior with canonical helper function
    _addPoseCanonical!( dfg, label, 0, prpo, genLabel=label, graphinit=graphinit, solvable=solvable,
                        srcType=varType, variableTags=variableTags, factorTags=factorTags,
                        postpose_cb=postpose_cb )
    #
  else
    @warn "$label already exists in the factor graph, no new variables added."
  end

  return dfg
end


"""
    $SIGNATURES

Build a basic factor graph in Pose2 with two `Pose2` and one landmark `Point2` variables,
along with `PriorPose2` on `:x0` and `Pose2Pose2` to `:x1`.  Also a `Pose2Point2BearingRange`
to landmark `:l1`.
"""
function generateCanonicalFG_TwoPoseOdo(;fg::AbstractDFG=initfg(),
                                        type::Type{<:Pose2}=Pose2,
                                        addlandmark::Bool=true,
                                        autoinit::Union{Bool,Nothing}=nothing,
                                        graphinit::Bool=true )
  #

  graphinit = if autoinit === nothing
    graphinit
  else
    @warn "autoinit is deprecated, use graphinit instead"
    autoinit
  end

  addVariable!(fg, :x0, Pose2)
  addVariable!(fg, :x1, Pose2)
  !addlandmark ? nothing : addVariable!(fg, :l1, Point2)

  addFactor!(fg, [:x0], PriorPose2(MvNormal([0;0;0.0],Matrix(Diagonal([1.0;1.0;0.01])))), graphinit=graphinit)
  addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal([10.0;0;0.0],Matrix(Diagonal([1.0;1.0;0.01])))), graphinit=graphinit)
  !addlandmark ? nothing : addFactor!(fg, [:x1;:l1], Pose2Point2BearingRange(Normal(0.0,0.01), Normal(20.0, 1.0)), graphinit=graphinit)

  return fg
end



#
