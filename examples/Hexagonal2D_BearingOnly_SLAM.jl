# import packages
using RoME, Distributions

# [OPTIONAL] add more julia processes to speed up inference
using Distributed
nprocs() < 3 ? addprocs(4-nprocs()) : nothing

@everywhere using RoME


# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), diagm( [0.1;0.1;0.1].^2 ))))

# Drive around in a hexagon a number of times
for i in 0:12
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], diagm( [0.1;0.1;0.1].^2 ) ))
    addFactor!(fg, [psym;nsym], pp )
end

# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
solveTree!(fg);


## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
using Cairo, RoMEPlotting
Gadfly.set_default_plot_size(35cm,20cm)


# For Juno/Jupyter style use
pl = plotSLAM2D(fg, contour=false)
# For scripting use-cases you can export the image
pl |> Gadfly.PDF("/tmp/test1.pdf", 20cm, 10cm)  # or PNG(...)


# Add a landmark l1
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])

bear1 = atan(10,20)
bear2 = atan(10,10) - pi/3


# Add landmarks with Bearing range measurements at x0, x6, x12, x18, x24, x30...
# Get variables
vars = ls(fg, r"x") |> sortDFG
for xIndex in 1:6:length(vars)
    @info("Creating factor between $(vars[xIndex]) and l1...")
    p2br = Pose2Point2Bearing(Normal(bear1,0.05))
    addFactor!(fg, [vars[xIndex]; :l1], p2br)
    p2br = Pose2Point2Bearing(Normal(bear2,0.05))
    addFactor!(fg, [vars[xIndex+1]; :l1], p2br)
end

# plotKDE(fg, [:x0;:x6;:x12], dims=[1;2])
# plotKDE(fg, [:x0;:x6;:x12], dims=[3])
# plotKDE(fg, [:x1;:x7;:x13], dims=[1;2])

# drawGraph(fg, show=true)

# Initialize :l1 numerical values but do not rerun solver
initAll!(fg)
pl = plotSLAM2D(fg, contour=false)
pl |> Gadfly.PDF("/tmp/test2.pdf", 20cm, 10cm)  # or PNG(...)

# solve, at time of writing compute time should less than 45s 
solveTree!(fg);

# redraw to see that additional information aids in refining the results
pl = plotSLAM2D(fg, contour=false)
pl |> Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm)  # or PNG(...)
