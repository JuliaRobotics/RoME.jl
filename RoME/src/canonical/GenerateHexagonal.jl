# generate canonical hexagonal graph


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

fg = generateGraph_Hexagonal()
drawGraph(fg, show=true)
```

Related

[`generateGraph_Circle`](@ref), [`generateGraph_Kaess`](@ref), [`generateGraph_TwoPoseOdo`](@ref), [`generateGraph_Boxes2D!`](@ref)
"""
function generateGraph_Hexagonal(;fg::AbstractDFG=initfg(),
                                  landmark::Bool=true,
                                  loopClosure::Bool=landmark,
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
  return generateGraph_Circle(6; graphinit, landmark, loopClosure, fg)
end
