# canonical factor graph examples useful for development and learning.

export loadCanonicalFG_Hexagonal, loadCanonicalFG_TwoPoseOdo
export warmUpSolverJIT

"""
    $SIGNATURES

Load a canonical factor graph: driving in a hexagonal circular pattern with one landmark.

Notes
- 7 Poses, :x0-:x6, Pose2,
- 1 Landmark, :l1, Point2,
- 6 Odometry, :x0x1f1, etc., Pose2Pose2 (Gaussian)
- 2 Sightings, :x0l1f1, :x6l1f1, RangeBearing (Gaussian)

Example
```julia
using RoME

fg = loadCanonicalFG_Hexagonal()
drawGraph(fg, show=true)
```
"""
function loadCanonicalFG_Hexagonal(;autoinit::Bool=true)
  # start with an empty factor graph object
  fg = initfg()
  # IIF.getSolverParams(fg).drawtree = true

  ## 1. Drive around in a hexagon
  # Add the first pose :x0
  println("STEP 1: Driving around a bit")
  addVariable!(fg, :x0, Pose2)
  # Add at a fixed location PriorPose2 to pin :x0 to a starting location
  prpo = PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I,3,3)))
  addFactor!(fg, [:x0], prpo, autoinit=autoinit)
  for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp , autoinit=autoinit)
  end

  # Add node linking initial pose with a bearing range measurement landmark
  addVariable!(fg, :l1, Point2, labels=[:LANDMARK;])
  p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
  addFactor!(fg, [:x0; :l1], p2br, autoinit=autoinit)

  # add loop closure sighting
  p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
  addFactor!(fg, [:x6; :l1], p2br, autoinit=autoinit)

  # return the new factor graph object
  return fg
end


"""
    $SIGNATURES

Build a basic factor graph in Pose2 with two `Pose2` and one landmark `Point2` variables,
along with `PriorPose2` on `:x0` and `Pose2Pose2` to `:x1`.  Also a `Pose2Point2BearingRange`
to landmark `:l1`.
"""
function loadCanonicalFG_TwoPoseOdo(;type::Type{Pose2}=Pose2,
                                    addlandmark::Bool=true,
                                    autoinit::Bool=true)
  #
  fg = initfg()

  addVariable!(fg, :x0, Pose2)
  addVariable!(fg, :x1, Pose2)
  !addlandmark ? nothing : addVariable!(fg, :l1, Point2)

  addFactor!(fg, [:x0], PriorPose2(MvNormal([0;0;0.0],Matrix(Diagonal([1.0;1.0;0.01])))), autoinit=autoinit)
  addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal([10.0;0;0.0],Matrix(Diagonal([1.0;1.0;0.01])))), autoinit=autoinit)
  !addlandmark ? nothing : addFactor!(fg, [:x1;:l1], Pose2Point2BearingRange(Normal(0.0,0.01), Normal(20.0, 1.0)), autoinit=autoinit)

  return fg
end


"""
    $SIGNATURES

Load and solve a canonical or user factor graph to warm up---precompile---several RoME/Caesar related functions.
"""
function warmUpSolverJIT(;drawtree::Bool=true)::Nothing
  #
  fg=loadCanonicalFG_Hexagonal()
  fcts = ls(fg, :x0)
  fcts = ls(fg)
  fcts = lsf(fg, :x0f1)
  fcts = lsf(fg)
  getSolverParams(fg).drawtree = drawtree
  tree, smt, hist = solveTree!(fg)
  nothing
end

#
