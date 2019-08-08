using RoME
using DistributedFactorGraphs

# start with an empty factor graph object
fg = LightGraphsDFG{SolverParams}(    params=SolverParams())

# Add the first pose :x0
x0 = addVariable!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
prior = addFactor!(fg, [:x0], PriorPose2( MvNormal([10; 10; 1.0/8.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) ) )

# Drive around in a hexagon in the cloud
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp )
end

# Alrighty! At this point, we should be able to solve locally...
# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
# Can do with graph too!
tree, smt, hist = solveTree!(fg)
