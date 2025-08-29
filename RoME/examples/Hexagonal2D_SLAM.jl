# add more julia processes
using Distributed
nprocs() < 3 ? addprocs(4-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using RoME
Distributed.@everywhere using RoME

# Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
using RoMEPlotting

## This part will be slow on first runs, also see pre-compiling docs page 
# warm up not-yet compiled code (optional step to speed up later parts)
# https://juliarobotics.org/Caesar.jl/latest/concepts/compile_binary/#Compiling-RoME.so
warmUpSolverJIT()


## build and solve the graph

fg = generateGraph_Hexagonal(landmark=false)


# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
tree = solveTree!(fg)

# also needs to compile first time
pl = plotSLAM2D(fg)
# For scripting use-cases you can export the image
pl |> Gadfly.PDF("/tmp/test1.pdf", 20cm, 10cm)  # or PNG(...)
# in REPL type ; to open shell, then view pdf, e.g. `evince /tmp/test1.pdf`

##

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, tags=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br )


# Initialize :l1 numerical values but do not rerun solver
initAll!(fg)
pl = plotSLAM2D(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test2.pdf", 20cm, 10cm),pl)  # or PNG(...)


## Add landmarks with Bearing range measurements

# add the loop closure
p2br2 = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br2 )

# and solve (which should be much faster now with all processes having necessary code compiled)
# again, also see: https://juliarobotics.org/Caesar.jl/latest/concepts/compile_binary/#Compiling-RoME.so
solveTree!(fg)

# redraw
pl = plotSLAM2D(fg, drawContour=false, drawEllipse=true)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)




#
