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


global _beehiveRecipe = Dict{Symbol,Symbol}(
  :l6   =>:l0,
  :l18  =>:l0,
  :l42  =>:l0,
  :l57  =>:l0,
  :l63  =>:l0,
  :l70  =>:l0,
  :l10  =>:l1,
  :l50  =>:l1,
  :l56  =>:l1,
  :l69  =>:l1,
  :l43  =>:l2,
  :l49  =>:l2,
  :l61  =>:l2,
  :l36  =>:l3,
  :l54  =>:l3,
  :l28  =>:l4,
  :l34  =>:l4,
  :l47  =>:l4,
  :l12  =>:l5,
  :l19  =>:l5,
  :l26  =>:l5,
  :l33  =>:l5,
  :l40  =>:l5,
  :l13  =>:l7,
  :l25  =>:l7,
  :l17  =>:l8,
  :l58  =>:l9,
  :l64  =>:l9,
  :l35  =>:l11,
  :l41  =>:l11,
  :l48  =>:l11,
  :l55  =>:l11,
  :l62  =>:l11,
  :l20  =>:l14,
  :l32  =>:l14,
  :l24  =>:l15,
  :l27  =>:l21,
  :l39  =>:l21,
  :l31  =>:l22,
  :l38  =>:l29,
  :l46  =>:l37,
  :l53  =>:l44,
  :l60  =>:l51,

  # pose offset legs
  :x41  =>:left,
  :x63  =>:left,
  :x78  =>:left,
)

# piecewise posecount construction
# 7,14,21,28,35,43,50



function _addLandmarkBeehive!(fg, 
                              lastPose::Symbol; 
                              refKey::Symbol=:simulated, 
                              solvable::Int=1,
                              graphinit::Bool=true  )
  #
  newFactor = RoME.Pose2Point2BearingRange(Normal(0,0.03), Normal(20,0.5))
  isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, lastPose, r"l\\d+", Point2, newFactor)

  # force isAlready until fixed _checkVariableByReference parametric solution at -pi Optim issue
  global _beehiveRecipe
  isAlready = false
  srcNumber = match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  genLabel = Symbol(:l, srcNumber)
  if haskey(_beehiveRecipe, genLabel) 
    isAlready = true
    genLabel = _beehiveRecipe[genLabel]
    simPPE = getPPE(fg, genLabel, refKey)
  end

  # maybe add new variable
  if !isAlready
    @info "New variable with simPPE" genLabel round.(simPPE.suggested,digits=2)
    newVar = addVariable!(fg, genLabel, Point2, solvable=solvable)
    addFactor!(fg, [lastPose; genLabel], newFactor, solvable=solvable, graphinit=graphinit)
    
    # also set :simulated PPE for similar future usage
    # newPPE = DFG.MeanMaxPPE(:simulated, simPPE, simPPE, simPPE)
    setPPE!(newVar, :simulated, typeof(simPPE), simPPE)   # TODO this API can be improved
  else
    @info "Adding simulated loop closure with perfect data association" lastPose genLabel
    addFactor!(fg, [lastPose; genLabel], newFactor, solvable=solvable, graphinit=graphinit)
  end

  #
  return genLabel
end

function _addPoseBeehive!(fg, 
                          prevLabel::Symbol,
                          posecount::Int,
                          factor::AbstractRelative;
                          refKey::Symbol=:simulated,
                          graphinit::Bool=false  )
  #
  # calculate and add the reference value
  isAlready, simPPE, = IIF._checkVariableByReference(fg, prevLabel, r"x\\d+", Pose2, factor, refKey=refKey)
  
  @show genLabel = Symbol("x", posecount)
  # add new pose variable
  v_n = addVariable!(fg, genLabel, Pose2 )
  addFactor!(fg, [prevLabel; genLabel], factor, graphinit=graphinit )

  # store simulated PPE for future use
  # ppe = DFG.MeanMaxPPE(refKey, simPPE, simPPE, simPPE)
  setPPE!(v_n, refKey, typeof(simPPE), simPPE)

  return genLabel
end


function _driveHex!(fgl::AbstractDFG,
                    posecount::Int;
                    poseCountTarget::Real=Inf,
                    graphinit::Bool=false,
                    gaugePrior::Symbol=:x0f1,
                    refKey::Symbol=:simulated,
                    addLandmarks::Bool=true,
                    landmarkSolvable::Int=1  )
  #

  # Drive around in a hexagon
  for i in (posecount):(posecount+5)
    if poseCountTarget <= posecount
      break
    end
    psym = Symbol("x$i")
    # nsym = Symbol("x$(i+1)")
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    posecount += 1
    nsym = _addPoseBeehive!(fgl, psym, posecount, pp, refKey=refKey, graphinit=graphinit)
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(fgl, nsym, refKey=refKey, solvable=landmarkSolvable, graphinit=false)
  end

  return posecount
end

function _offsetHexLeg( fgl::AbstractDFG,
                        posecount::Int; 
                        poseCountTarget::Real=Inf,
                        direction::Symbol=:right,
                        graphinit::Bool=false,
                        psym = Symbol("x$(posecount)"),
                        # nsym = Symbol("x$(posecount+1)"),
                        refKey::Symbol=:simulated,
                        guagePrior::Symbol=:x0f1,
                        addLandmarks::Bool=true,
                        landmarkSolvable::Int=1   )
  #
  #skip out early if pose count target has been reached
  if poseCountTarget <= posecount
    return posecount
  end
  dirsign = if direction == :right
    -1
  elseif direction == :left
    +1
  else
    error("unknown direction symbol $direction")
  end
  pp = Pose2Pose2(MvNormal([10.0;0;dirsign*pi/3], diagm([0.1;0.1;0.1].^2)))
  
  posecount += 1
  nsym = _addPoseBeehive!(fgl, psym, posecount, pp, refKey=refKey, graphinit=graphinit)

  # add a new landmark (if not yet present)
  !addLandmarks ? nothing : _addLandmarkBeehive!(fgl, nsym, refKey=refKey, solvable=landmarkSolvable, graphinit=false)

  return posecount
end


function generateCanonicalFG_Beehive!(poseCountTarget::Int=36;
                                      fg::AbstractDFG = initfg(),
                                      direction::Symbol = :right,
                                      graphinit::Bool = false,
                                      refKey::Symbol=:simulated,
                                      addLandmarks::Bool=true,
                                      landmarkSolvable::Int=0,
                                      useMsgLikelihoods::Bool=false     )
  #
  global _beehiveRecipe

  # does anything exist in the graph yet
  posecount = if :x0 in ls(fg)
    # what is the last pose
    lastPose = (ls(fg, r"x\d+") |> sortDFG)[end]
    # get latest posecount number
    match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  else
    # initial zero pose
    generateCanonicalFG_ZeroPose2(fg=fg, graphinit=graphinit) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue
    getSolverParams(fg).useMsgLikelihoods = useMsgLikelihoods    

    # reference ppe on :x0
    refVal = zeros(3)
    ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
    setPPE!(fg[:x0], refKey, DFG.MeanMaxPPE, ppe)
    
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(fg, :x0, refKey=refKey, solvable=landmarkSolvable, graphinit=false)

    # staring posecount (i.e. :x0)
    0
  end

  # keep adding poses until the target number is reached
  while posecount < poseCountTarget
    posecount = _driveHex!(fg, posecount, graphinit=graphinit, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget)
    # drive the offset legs
    lastPose = Symbol(:x, posecount)
    if haskey(_beehiveRecipe, lastPose)
      posecount = _offsetHexLeg(fg, posecount, direction=_beehiveRecipe[lastPose], graphinit=graphinit, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget)    
    end
    posecount = _offsetHexLeg(fg, posecount, direction=direction, graphinit=graphinit, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget)
  end

  # FIXME force everything solvable for now
  setSolvable!.(fg, ls(fg), 1)
  setSolvable!.(fg, lsf(fg), 1)

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
