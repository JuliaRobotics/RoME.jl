# canonical factor graph examples useful for development and learning.

export generateCanonicalFG_ZeroPose2, generateCanonicalFG_Hexagonal, generateCanonicalFG_TwoPoseOdo, generateCanonicalFG_Circle
export warmUpSolverJIT

"""
    $SIGNATURES

Generate a canonical factor graph with a Pose2 `:x0` and MvNormal with covariance `P0`
"""
function generateCanonicalFG_ZeroPose2(;fg::AbstractDFG=initfg(),
                                       graphinit::Bool=true,
                                       Σ0::AbstractMatrix{<:Real}= 0.01*Matrix{Float64}(LinearAlgebra.I,3,3),
                                       μ0::AbstractVector{<:Real}=zeros(3),
                                       label::Symbol=:x0 )
  #
  if !exists(fg, label)
    addVariable!(fg, label, Pose2)
    prpo = PriorPose2(MvNormal(μ0, Σ0))
    addFactor!(fg, [label], prpo, graphinit=graphinit)
  end
  return fg
end

"""
    $SIGNATURES

Generate a canonical factor graph: driving in a circular pattern with one landmark.

Notes
- Poses, :x0, :x1,... Pose2,
- Odometry, :x0x1f1, etc., Pose2Pose2 (Gaussian)
- OPTIONAL: 1 Landmark, :l1, Point2,
- 2 Sightings, :x0l1f1, :x6l1f1, RangeBearing (Gaussian)

Example
```julia
using RoME

fg = generateCanonicalFG_Hexagonal()
drawGraph(fg, show=true)
```

Related

generateCanonicalFG_Circle, generateCanonicalFG_Kaess, generateCanonicalFG_TwoPoseOdo
"""
function generateCanonicalFG_Circle(poses::Int=6;
                                    fg::AbstractDFG=initfg(),
                                    offsetPoses::Int=maximum([length(ls(fg, r"x\d"))-1;0]),
                                    autoinit::Union{Bool, Nothing}=nothing,
                                    graphinit::Bool=true,
                                    landmark::Bool=true,
                                    loopClosure::Bool=true,
                                    stopEarly::Int=9999999,
                                    biasTurn::Real=0.0,
                                    kappaOdo::Real=1.0,
                                    cyclePoses::Int=poses )
  # assume empty factor graph object fg
  @assert offsetPoses < poses "`offsetPoses` must be smaller than total number of `poses`"
  # IIF.getSolverParams(fg).drawtree = true

  graphinit = if autoinit === nothing
    graphinit
  else
    @warn "autoinit is deprecated, use graphinit instead"
    autoinit
  end

  ## 1. Drive around in a hexagon
  # Add the first pose :x0
  if !exists(fg, :x0)
    addVariable!(fg, :x0, Pose2)
    prpo = PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I,3,3)))
    addFactor!(fg, [:x0], prpo, graphinit=graphinit)
  end

  # println("STEP 1: Driving around a bit")
  # Add at a fixed location PriorPose2 to pin :x0 to a starting location
  for i in offsetPoses:poses-1
    if stopEarly <= i
      break
    end
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;2pi/(cyclePoses)+biasTurn], Matrix(Diagonal((kappaOdo*[0.1;0.1;0.1]).^2))))
    addFactor!(fg, [psym;nsym], pp , graphinit=graphinit)
  end

  if !landmark
    return fg
  end
  if !exists(fg, :l1)
    # Add node linking initial pose with a bearing range measurement landmark
    addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
    p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
    addFactor!(fg, [:x0; :l1], p2br, graphinit=graphinit)
  end

  if !loopClosure || !exists(fg, Symbol("x$poses"))
    return fg
  end
  # add loop closure sighting
  p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
  addFactor!(fg, [Symbol("x$poses"); :l1], p2br, graphinit=graphinit)

  # return the new factor graph object
  return fg
end



"""
    $SIGNATURES

Generate a canonical factor graph: driving in a hexagonal circular pattern with one landmark.

Notes
- 7 Poses, :x0-:x6, Pose2,
- 1 Landmark, :l1, Point2,
- 6 Odometry, :x0x1f1, etc., Pose2Pose2 (Gaussian)
- 2 Sightings, :x0l1f1, :x6l1f1, RangeBearing (Gaussian)

Example
```julia
using RoME

fg = generateCanonicalFG_Hexagonal()
drawGraph(fg, show=true)
```

Related

generateCanonicalFG_Circle, generateCanonicalFG_Kaess, generateCanonicalFG_TwoPoseOdo
"""
function generateCanonicalFG_Hexagonal(;fg::AbstractDFG=initfg(),
                                        N::Int=100,
                                        autoinit::Union{Bool, Nothing}=nothing,
                                        graphinit::Bool=true  )
  #
  getSolverParams(fg).N = N
  graphinit = if autoinit === nothing
    graphinit
  else
    @warn "autoinit is deprecated, use graphinit instead"
    autoinit
  end
  return generateCanonicalFG_Circle(6, graphinit=graphinit, landmark=true, loopClosure=true; fg=fg)
end


"""
    $SIGNATURES

If you know a variable is `::Type{<:Pose2}` but want to find its default prior `::Type{<:PriorPose2}`.

Assumptions
- The prior type will be defined in the same module as the variable type.
- Not exported per default, but can be used with knowledge of the caveats.

Example
```julia
using RoME
@assert RoME.PriorPose2 == DFG._getPriorType(Pose2)
```

DevNotes
- TODO move to DFG instead
"""
_getPriorType(_type::InferenceVariable) = getfield(_type.name.module, Symbol(:Prior, _type.name.name))


"""
    $SIGNATURES

Check if a variable might already be located at the test location, by means of a (default) `refKey=:simulated` PPE stored in the existing variables.

Notes
- Checks, using provided `factor` from `srcLabel` in `fg` to an assumed `dest` variable whcih may or may not yet exist.
- This function was written to aid in building simulation code, 
  - it's use in real world usage may have unexpected behaviour -- hence not exported.
- Return `::Tuple{Bool, Vector{Float64}, Symbol}`, eg. already exists `(true, [refVal], :l17)`, or if a refernce variable does not yet `(false, [refVal], :l28)`.
  - Vector contains the PPE reference location of the new variable as calculated.
- Auto `destPrefix` is trying to parse `destRegex` labels like `l\\d+` or `tag\\d+`, won't work with weirder labels e.g. `:l_4_23`.
  - User can overcome weird names by self defining `destPrefix` and `srcNumber`.
  - User can also ignore and replace the generated new label `Symbol(destPrefix, srcNumber)`.
- This function does not add new variables or factors to `fg`, user must do that themselves after.
  - Useful to use in combination with `setPPE!` on new variable.
- At time of writing `accumulateFactorMeans` could only incorporate priors or binary relative factors.
  - internal info, see [`solveBinaryFactorParameteric`](@ref),
  - This means at time of writing `factor` must be a binary factor.
- Tip, if simulations are inducing odometry bias, think of using two factors from caller (e.g. simPerfect and simBias).

Example
```julia
# fg has :x5 and :l2 and PPEs :simulated exists in all variables
# user wants to add a factor from :x5 to potential new :l5, but maybe a (simulated) variable, say :l2, is already there.

newFactor = RoME.Pose2Point2BearingRange(Normal(), Normal(20,0.5))
isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, :x5, r"l\\d+", Point2, newFactor)

# maybe add new variable
if !isAlready
  @info "New variable with simPPE" genLabel simPPE 
  newVar = addVariable!(fg, genLabel, Point2)
  addFactor!(fg, [:x5; genLabel], newFactor)

  # also set :simulated PPE for similar future usage
  newPPE = DFG.MeanMaxPPE(:simulated, simPPE, simPPE, simPPE)
  setPPE!(newVar, :simulated, typeof(newPPE), newPPE)   # TODO this API can be improved
else
  @info "Adding simulated loop closure with perfect data association" :x5 genLabel
  addFactor!(fg, [:x5; genLabel], newFactor)
end

# the point is that only the (0,20) values in newFactor are needed, all calculations are abstracted away.
```


Related

[`RoME.generateCanonicalFG_Beehive!`](@ref), [`accumulateFactorMeans`](@ref), [`getPPE`](@ref)
"""
function _checkVariableByReference( fg::AbstractDFG,
                                    srcLabel::Symbol,            # = :x5
                                    destRegex::Regex,            # = r"l\d+"
                                    destType::InferenceVariable, # = Point2
                                    factor::AbstractRelative;    # = Pose2Poin2BearingRange(...)
                                    srcType::InferenceVariable = getVariableType(fg, srcLabel) |> typeof,
                                    refKey::Symbol=:simulated,
                                    prior = _getPriorType(srcType)( MvNormal(getPPE(fg[srcLabel], refKey).suggested, diagm(ones(getDimension(srcType)))) ),
                                    atol::Real = 1e-3,
                                    destPrefix::Symbol = match(r"[a-zA-Z]+", destRegex.pattern) |> Symbol,
                                    srcNumber = match(r"\d+", string(srcLabel)).match |> x->parse(Int,x)  )
  #
  
  # calculate and add the reference value
  lPp = 
  tfg = initfg()
  addVariable!(tfg, :x0, srcType )
  addFactor!(tfg, [:x0], prior )
  addVariable!(tfg, :l0, destType )
  addFactor!( tfg, [:x0; :l0], factor, graphinit=false )
  
  # calculate where the landmark reference position is
  refVal = accumulateFactorMeans(tfg, [:x0f1; :x0l0f1])
  ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
  # setPPE!(v_n, refKey, DFG.MeanMaxPPE, ppe)

  # now check if we already have a landmark at this location
  varLms = ls(fg, destRegex)
  ppeLms = getPPE.(getVariable.(fg, varLms), refKey) .|> x->x.suggested
  errmask = ppeLms .|> x -> norm(x - ppe) < atol
  already = any(errmask)

  alrLm = varLms[findfirst(errmask)]
  @assert sum(errmask) == 1 "Besides $alrLm, there should be only one landmark at $ppe"
  if already
    # does exist, ppe, variableLabel
    return true, ppe, alrLm
  end
  
  # Nope does not exist, ppe, generated new variable label only
  return false, ppe, Symbol(destPrefix, srcNumber)
end

function _addLandmarkBeehive!(fg, lastPose::Symbol)
  #
  newFactor = RoME.Pose2Point2BearingRange(Normal(0,0.03), Normal(20,0.5))
  isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, :x5, r"l\\d+", Point2, newFactor)

  # maybe add new variable
  if !isAlready
    @info "New variable with simPPE" genLabel simPPE 
    newVar = addVariable!(fg, genLabel, Point2)
    addFactor!(fg, [lastPose; genLabel], newFactor)

    # also set :simulated PPE for similar future usage
    newPPE = DFG.MeanMaxPPE(:simulated, simPPE, simPPE, simPPE)
    setPPE!(newVar, :simulated, typeof(newPPE), newPPE)   # TODO this API can be improved
  else
    @info "Adding simulated loop closure with perfect data association" lastPose genLabel
    addFactor!(fg, [lastPose; genLabel], newFactor)
  end

  #
  return genLabel
end


function _driveHex!(fgl::AbstractDFG,
                    posecount::Int;
                    graphinit::Bool=false,
                    gaugePrior::Symbol=:x0f1,
                    refKey::Symbol=:simulated )
  #

  # Drive around in a hexagon
  for i in (posecount):(posecount+5)
    psym = Symbol("x$i")
    posecount += 1
    nsym = Symbol("x$(i+1)")
    v_n = addVariable!(fgl, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fgl, [psym;nsym], pp, graphinit=graphinit )

    # add a new landmark
    _addLandmarkBeehive!(fgl, nsym)

    # # calculate and add the reference value
    # pth = findShortestPathDijkstra(fgl, gaugePrior, nsym, regexVariables=r"x\d+")
    # pth = pth[isFactor.(fgl, pth)]
    # refVal = accumulateFactorMeans(fgl, pth)
    # ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
    # setPPE!(v_n, refKey, DFG.MeanMaxPPE, ppe)
  end

  return posecount
end

function _offsetHexLeg( fgl::AbstractDFG,
                        posecount::Int; 
                        direction::Symbol=:right,
                        graphinit::Bool=false,
                        psym = Symbol("x$(posecount)"),
                        nsym = Symbol("x$(posecount+1)"),
                        refKey::Symbol=:simulated,
                        guagePrior::Symbol=:x0f1   )
  #
  posecount += 1
  v_n = addVariable!(fgl, nsym, Pose2 )
  pp = nothing
  if direction == :right
      pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  elseif direction == :left
      pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  end
  addFactor!(fgl, [psym; nsym], pp, graphinit=graphinit )

  # calculate and add the reference value
  pth = findShortestPathDijkstra(fgl, gaugePrior, nsym, regexVariables=r"x\d+")
  pth = pth[isFactor.(fgl, pth)]
  refVal = accumulateFactorMeans(fgl, pth)
  ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
  setPPE!(v_n, refKey, DFG.MeanMaxPPE, ppe)

  return posecount
end


function generateCanonicalFG_Beehive!(cells::Int=2;
                                      fg::AbstractDFG = initfg(),
                                      direction::Symbol = :right,
                                      graphinit::Bool = false,
                                      refKey::Symbol=:simulated    )
  #
  
  # does anything exist in the graph yet
  posecount = if :x0 in ls(fg)
    # what is the last pose
    lastPose = (ls(fg, r"x\d+") |> sortDFG)[end]
    # get latest posecount number
    match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  else
    # initial zero pose
    generateCanonicalFG_ZeroPose2(fg=fg, graphinit=graphinit)
    
    # Add landmarks with Bearing range measurements
    addVariable!(fg, :l0, Point2, tags=[:LANDMARK;])
    p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
    addFactor!(fg, [:x0; :l0], p2br, graphinit=false )
    
    # reference ppe on :x0 and :l0
    refVal = zeros(3)
    ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
    setPPE!(fg[:x0], refKey, DFG.MeanMaxPPE, ppe)
    refVal = accumulateFactorMeans(fgl, [:x0f1; :x0l0f1])
    ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
    setPPE!(fg[:l0], refKey, DFG.MeanMaxPPE, ppe)
    # staring posecount (i.e. :x0)
    0
  end
  
  posecount = _driveHex!(fg, posecount)
  # drive the offset leg
  posecount = _offsetHexLeg(fg, posecount, direction=direction, graphinit=graphinit)
  # drive the next hex
  
  return fg
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



"""
    $SIGNATURES

Load and solve a canonical or user factor graph to warm up---precompile---several RoME/Caesar related functions.
"""
function warmUpSolverJIT(;fg::AbstractDFG=generateCanonicalFG_Hexagonal(),
                          drawtree::Bool=true )::Nothing
  #

  fcts = ls(fg, :x0)
  fcts = ls(fg)
  fcts = lsf(fg, :x0f1)
  fcts = lsf(fg)
  getSolverParams(fg).drawtree = drawtree
  tree, smt, hist = solveTree!(fg)
  nothing
end

#
