# generate canonical circular graph



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

DevNotes
- TODO refactor to use new `_calcHelix_T`.

Related

[`generateCanonicalFG_Circle`](@ref), [`generateCanonicalFG_Kaess`](@ref), [`generateCanonicalFG_TwoPoseOdo`](@ref)
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
    # TODO rework the addition of noise alongside simulated ground truth, maybe `noise_t_cb`.
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