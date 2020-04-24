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

  graphinit = if autoinit == nothing
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

  println("STEP 1: Driving around a bit")
  # Add at a fixed location PriorPose2 to pin :x0 to a starting location
  for i in offsetPoses:poses-1
    if stopEarly <= i
      break
    end
    @show psym = Symbol("x$i")
    @show nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;2pi/(cyclePoses)+biasTurn], Matrix(Diagonal((kappaOdo*[0.1;0.1;0.1]).^2))))
    addFactor!(fg, [psym;nsym], pp , graphinit=graphinit)
  end

  if !landmark
    return fg
  end
  if !exists(fg, :l1)
    # Add node linking initial pose with a bearing range measurement landmark
    addVariable!(fg, :l1, Point2, labels=[:LANDMARK;])
    p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
    addFactor!(fg, [:x0; :l1], p2br, graphinit=graphinit)
  end

  if !loopClosure || !exists(fg, :x*"$poses")
    return fg
  end
  # add loop closure sighting
  p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
  addFactor!(fg, [:x*"$poses"; :l1], p2br, graphinit=graphinit)

  # return the new factor graph object
  return fg
end


function driveHex(fgl, posecount::Int; steps::Int=5)
    # Drive around in a hexagon
    for i in (posecount-1):(posecount-1+steps)
        psym = Symbol("x$i")
        posecount += 1
        nsym = Symbol("x$(i+1)")
        addVariable!(fgl, nsym, Pose2)
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
        addFactor!(fgl, [psym;nsym], pp, autoinit=false )
    end

    return posecount
end


function offsetHexLeg(dfg::G, posecount::Int; direction=:right) where G <: AbstractDFG
    psym = Symbol("x$(posecount-1)")
    nsym = Symbol("x$(posecount)")
    posecount += 1
    addVariable!(dfg, nsym, Pose2)
    pp = nothing
    if direction == :right
        pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    elseif direction == :left
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    end
    addFactor!(dfg, [psym; nsym], pp, autoinit=false )
    return posecount
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
                                       autoinit::Union{Bool, Nothing}=nothing,
                                       graphinit::Bool=true  )
  #
  graphinit = if autoinit == nothing
    graphinit
  else
    @warn "autoinit is deprecated, use graphinit instead"
    autoinit
  end
  return generateCanonicalFG_Circle(6, graphinit=graphinit, landmark=true, loopClosure=true)
end


"""
    $SIGNATURES

Build a basic factor graph in Pose2 with two `Pose2` and one landmark `Point2` variables,
along with `PriorPose2` on `:x0` and `Pose2Pose2` to `:x1`.  Also a `Pose2Point2BearingRange`
to landmark `:l1`.
"""
function generateCanonicalFG_TwoPoseOdo(;fg::AbstractDFG=initfg(),
                                        type::Type{Pose2}=Pose2,
                                        addlandmark::Bool=true,
                                        autoinit::Union{Bool,Nothing}=nothing,
                                        graphinit::Bool=true )
  #

  graphinit = if autoinit == nothing
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
