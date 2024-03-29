# canonical factor graph examples useful for development and learning.


"""
    $SIGNATURES

Lower level utility function for simulated canonical graph generation, where a pose variable and factor is added using this function.

Notes
- Also adds the `refKey=:simulated` (default) PPE as in graph record of the simulated values. 
- Use callback `postpose_cb(g::AbstractDFG,lastpose::Symbol)` to call user operations after each pose step.

Related

[`generateGraph_ZeroPose`](@ref), [`generateGraph_Helix2D!`](@ref)
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
                            inflation::Real=getSolverParams(fg).inflation,
                            variableTags::AbstractVector{Symbol}=Symbol[],
                            factorTags::AbstractVector{Symbol}=Symbol[],
                            doRef::Bool=true,
                            refKey::Symbol = :simulated,
                            overridePPE = doRef ? nothing : zeros(getDimension(poseType)),
                            postpose_cb::Function=(fg_,latestpose)->()  )
  #
  # calculate and add the reference value
  isAlready, simPPE, = IIF._checkVariableByReference(fg, prevLabel, poseRegex, poseType, factor, refKey=refKey, srcType=srcType, doRef=doRef, overridePPE=overridePPE)

  # dispatch on prior or binary factor
  _getlabels(fact::AbstractPrior) = [genLabel;]
  _getlabels(fact::AbstractRelative) = [prevLabel; genLabel]
  
  # add new pose variable
  v_n = addVariable!(fg, genLabel, poseType; solvable, tags=variableTags )
  addFactor!(fg, _getlabels(factor), factor; graphinit, solvable, tags=factorTags, inflation )

  # store simulated PPE for future use
  # ppe = DFG.MeanMaxPPE(refKey, simPPE, simPPE, simPPE)
  doRef ? setPPE!(v_n, refKey, typeof(simPPE), simPPE) : nothing

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
function generateGraph_ZeroPose(; varType::Type{<:InferenceVariable}=Pose2,
                                  graphinit = nothing,
                                  solverParams::SolverParams=SolverParams(),
                                  dfg::AbstractDFG = LocalDFG{SolverParams}(;solverParams),  
                                  doRef::Bool=true,
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
  (graphinit isa Nothing) ? nothing : @error("generateGraph_ZeroPose keyword graphinit obsolete, use solverParams(graphinit=..) instead.")

  # only add the first variable if none others exist
  if !exists(dfg, label)
    # generate a default prior
    prpo = priorType(priorArgs...)
    # add the variable and prior with canonical helper function
    _addPoseCanonical!( dfg, label, 0, prpo; genLabel=label, solvable,
                        srcType=varType, variableTags, factorTags, graphinit = solverParams.graphinit,
                        doRef, postpose_cb )
    #
  else
    @warn "$label already exists in the factor graph, no new variables added."
  end

  return dfg
end


"""
    $SIGNATURES

Build a chain of variables with binary factors.

Notes:
- `fctData[1]` aligns with the first variable (likely `:x0` for default empty `dfg`), and 
  - `:x0` has no previous variable/factor data.
  
DevNotes
- TODO Consolidate with [`IIF._buildGraphByFactorAndTypes!`](@ref) and wrap in RoME?
"""
function buildGraphChain!(  fctData::AbstractVector = [MvNormal([10;0;0.],diagm(0.1.*ones(3))) for _ in 1:3],
                            fctType::Type{<:AbstractRelative} = Pose2Pose2,
                            preFct_args_cb::Function = (fg_,data)->(data.currData,);
                            stopAfter::Integer=2^(Sys.WORD_SIZE-1)-1,
                            varType::IIF.InstanceType{InferenceVariable} = Pose2,
                            solverParams::SolverParams = SolverParams(),
                            solvable::Integer = 1,
                            fctKwargs::NamedTuple = (;),
                            varTags::AbstractVector{Symbol} = Symbol[],
                            fctTags::AbstractVector{Symbol} = Symbol[],
                            postpose_cb::Function = (fg_, lastpose) -> nothing,
                            doRef::Bool=true,
                            dfg::AbstractDFG = generateGraph_ZeroPose(; varType, solverParams, solvable, doRef, postpose_cb),
                            inflation_fct::Real = getSolverParams(dfg).inflation,
                            varRegex::Regex = r"x\d+",
                            varLast::Symbol = sortDFG(ls(dfg, varRegex))[end],
                            varCount::Integer = match(r"\d+", string(varLast)).match |> x->parse(Int,x),
                            varPrefix::Symbol = match(r"[a-zA-Z_]+", varRegex.pattern).match |> Symbol  )
  #

  for i in 1:length(fctData)
    i <= stopAfter ? nothing : break
    varCount += 1
    varCurr = Symbol(varPrefix, varCount)
    # prep small container for data
    prev_ = 1 <= i-1 ? fctData[i-1] : nothing
    next_ = i+1 <= length(fctData) ? fctData[i+1] : nothing
    dat_ = (;varCurr=varCurr, varLast=varLast, prevData=prev_, currData=fctData[i], nextData=next_)
    _addPoseCanonical!( dfg, 
                        varLast, 
                        varCount,
                        fctType(preFct_args_cb(dfg,dat_)...; fctKwargs...);
                        doRef=doRef,
                        poseRegex=varRegex,
                        genLabel=varCurr,
                        srcType=varType,
                        solvable=solvable,
                        inflation=inflation_fct,
                        variableTags=varTags,
                        factorTags=fctTags,
                        postpose_cb=postpose_cb )
    #
    varLast = varCurr
  end
  
  dfg
end



"""
    $SIGNATURES

Build a basic factor graph in Pose2 with two `Pose2` and one landmark `Point2` variables,
along with `PriorPose2` on `:x0` and `Pose2Pose2` to `:x1`.  Also a `Pose2Point2BearingRange`
to landmark `:l1`.

DevNotes
- TODO standardize to latest and greatest generator api pattern

See also: [`buildGraphChain!`](@ref)
"""
function generateGraph_TwoPoseOdo(; solverParams::SolverParams = SolverParams(),
                                    varType::Type{<:Pose2}=Pose2,
                                    solvable::Integer = 1,
                                    doRef::Bool=true,
                                    postpose_cb::Function = (fg_, lastpose) -> nothing,
                                    dfg::AbstractDFG=generateGraph_ZeroPose(; varType, solverParams, solvable, doRef, postpose_cb),
                                    addlandmark::Bool=true,
                                    fg=nothing,
                                    graphinit=nothing )
  #
  (graphinit isa Nothing) ? nothing : @error("generateGraph_TwoPoseOdo keyword graphinit obsolete, use solverParams=SolverParams(graphinit=..) instead.")
  (fg isa Nothing) ? nothing : @error("generateGraph_TwoPoseOdo keyword fg obsolete, use dfg= instead.")

  buildGraphChain!( [MvNormal([10.0;0;0.0],diagm([1.0;1.0;0.01]));];
                    dfg,
                    varType )
  
  if addlandmark 
    addVariable!(dfg, :l1, Point2)
    addFactor!( dfg, [:x1;:l1], Pose2Point2BearingRange(Normal(0.0,0.01), Normal(20.0, 1.0)), graphinit=solverParams.graphinit )
  end
  
  # return the built graph
  return dfg
end

"""
    $SIGNATURES

Simulates IMU measurements from body frame rates and desired world frame accelerations.
"""
function generateField_InertialMeasurement(;
  dt = 0.01,
  N = 401,
  rate = [0.0, 0.0, pi/2],
  w_R_b = [1. 0 0; 0 1 0; 0 0 1],
  gravity = [0.0, 0, 0],
  accel0 = [0.0, 0, 0] + gravity,
  b_a = SA[0.0, 0, 0], # [0.0, pi/2*10, 0],
  σ_a = 0.0, # 1e-4, #0.16e-3*9.81  # noise density m/s²/√Hz
  σ_ω = 0.0, # deg2rad(0.0001),  # noise density rad/√Hz
)
  tspan = (0.0, dt*(N-1))
  
  gn = norm(σ_ω) < 1e-14 ? ()->[0, 0, 0] : ()->rand(MvNormal(diagm(ones(3)*σ_ω^2 * 1/dt)))
  an = norm(σ_a) < 1e-14 ? ()->[0, 0, 0] : ()->rand(MvNormal(diagm(ones(3)*σ_a^2 * 1/dt)))

  Σy  = diagm([ones(3)*σ_a^2; ones(3)*σ_ω^2])

  gyros = [rate + gn() for _ = 1:N]
  
  accels = Vector{typeof(accel0)}()
  push!(accels, deepcopy(accel0) + an())
  # accels = [deepcopy(accel0) + an()]
  M = SpecialOrthogonal(3)

  # b_a = [0.1, 0, 0]
  for g in gyros[1:end-1]
    X = hat(M, Identity(M), g)
    exp!(M, w_R_b, w_R_b, X*dt)
    push!(accels, (b_a .+ an()) + w_R_b' * accel0)
  end

  (;tspan,gyros,accels,Σy)
end

generateField_InertialMeasurement_RateZ(;
  dt = 0.01,
  N = 401,
  rate = [0.0, 0.0, pi/2],
  kw...
) = generateField_InertialMeasurement(;dt,N,rate,kw...)


generateField_InertialMeasurement_noise(;
  dt = 0.1,
  N = 11,
  rate = [0, 0, 0.001],
  gravity = SA[0, 0, 9.81],
  accel0 = [0, 0, -1.0] + gravity,
  σ_a = 1e-4,             # 0.16e-3*9.81  # noise density m/s²/√Hz
  σ_ω = deg2rad(0.0001),  # noise density rad/√Hz
) = generateField_InertialMeasurement(;
  dt, 
  N, 
  rate, 
  gravity, 
  accel0, 
  σ_a, 
  σ_ω
)





#
