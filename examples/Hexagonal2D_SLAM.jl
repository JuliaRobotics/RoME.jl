# add more julia processes
nprocs() < 3 ? addprocs(4-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using RoME, Distributions


# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addNode!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*eye(3))))

# Drive around in a hexagon
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addNode!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], diagm([0.1;0.1;0.1].^2)))
  addFactor!(fg, [psym;nsym], pp )
end

# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)
# batchSolve!(fg) # coming soon



## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
using RoMEPlotting

# For Juno/Jupyter style use
pl = drawPoses(fg)
# For scripting use-cases you can export the image
Gadfly.draw(Gadfly.PDF("/tmp/test1.pdf", 20cm, 10cm),pl)  # or PNG(...)

# Add landmarks with Bearing range measurements
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br)

# Initialize :l1 numerical values but do not rerun solver
ensureAllInitialized!(fg)
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test2.pdf", 20cm, 10cm),pl)  # or PNG(...)


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br2)

# solve
tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)

# redraw
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)


# using RoMEPlotting, Gadfly
#
#
#
#
#
# pl = plotKDE(fg, [:x0; :x1; :x2; :x3; :x4; :x5; :x6]);
#
# Gadfly.draw(PDF("tmpX0123456.pdf", 15cm, 15cm), pl)
#
# @async run(`evince tmpX0123456.pdf`)
#
#
#
# # pl = drawPoses(fg)
# pl = drawPosesLandms(fg)
# Gadfly.draw(PDF("tmpPosesFg.pdf", 15cm, 15cm), pl)
# @async run(`evince tmpPosesFg.pdf`)
#
#
#
# tree = wipeBuildNewTree!(fg)
#
#
# @async Graphs.plot(tree.bt)
#
#
# @time inferOverTree!(fg, tree)
#
#
# # Graphs.plot(tree.bt)
#
# @async Graphs.plot(fg.g)
#
#
# # These functions need more work
#
# using KernelDensityEstimate
#
# function plot(::Type{Pose2}, p::Vector{BallTreeDensity})
#
#   xy = plotKDE(p,dims=[1;2])
#   th = plotKDE(p,dims=[3])
#
#   hstack(xy, th)
# end
#
# vert0 = getVert(fg, :x0)
# vert1 = getVert(fg, :x1)
# pl = plot(typeof(getData(vert).softtype), [getKDE(vert0); getKDE(vert1)])
#
# Gadfly.draw(PDF("tmpX01.pdf", 15cm, 15cm), pl)
# @async run(`evince tmpX01.pdf`)


#
