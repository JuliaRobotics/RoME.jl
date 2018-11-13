# add more julia processes
using Distributed
nprocs() < 3 ? addprocs(4-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using IncrementalInference
using RoME, Distributions
## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
using RoMEPlotting


# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addNode!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], IIF.Prior( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2)) )))

# Drive around in a hexagon
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addNode!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  addFactor!(fg, [psym;nsym], pp )
end

# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
batchSolve!(fg)

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
batchSolve!(fg)

# redraw
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)
