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

# Drive around in a hexagon a number of times
for i in 0:12
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addNode!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], diagm([0.1;0.1;0.1].^2)))
    addFactor!(fg, [psym;nsym], pp )
end

# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
batchSolve!(fg)

## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
using RoMEPlotting
using Compose

# For Juno/Jupyter style use
pl = drawPoses(fg)
# For scripting use-cases you can export the image
Gadfly.push_theme(:default) # light background, where Juno uses dark background
Gadfly.draw(Gadfly.PDF("/tmp/test1.pdf", 20cm, 10cm),pl)  # or PNG(...)

# Add a landmark l1
addNode!(fg, :l1, Point2, labels=["LANDMARK"])

# Add landmarks with Bearing range measurements at x0, x6, x12, x18, x24, x30...
vars = ls(fg)[1] # Get variables
for xIndex in 1:6:length(vars)
    info("Creating factor between $(vars[xIndex]) and l1...")
    p2br = Pose2Point2Bearing(Normal(0,0.05))
    addFactor!(fg, [vars[xIndex]; :l1], p2br)
end

# Initialize :l1 numerical values but do not rerun solver
ensureAllInitialized!(fg)
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test2.pdf", 20cm, 10cm),pl)  # or PNG(...)

# solve
batchSolve!(fg)

# redraw to see that additional information aids in refining the results
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)
