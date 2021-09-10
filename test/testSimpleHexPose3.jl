using RoME
using Test
using DistributedFactorGraphs
# @testset "Test hexagonal for pose3" begin

# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
x0 = addVariable!(fg, :x0, Pose3)

# Add at a fixed location PriorPose3 
prior = addFactor!(fg, [:x0], PriorPose3( MvNormal([0., 0, 0, 0, 0, 0], [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]))) 

# Drive around in a hexagon horizontally
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose3)
    pp = Pose3Pose3(MvNormal([10.0, 0, 0, 0, 0, pi/3], [0.5, 0.5, 0.5, 0.05, 0.05, 0.05]))
    addFactor!(fg, [psym;nsym], pp )
end

ensureAllInitialized!(fg)
getPPESuggested.(fg, sortDFG(ls(fg)))

# tree = solveTree!(fg)


# end

# @testset "Test hexagonal for pose3" begin

# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
x0 = addVariable!(fg, :x0, Pose3)

# Add at a fixed location PriorPose3 
# prior = addFactor!(fg, [:x0], PriorPose3( MvNormal([2., 2, 2, 0, 0, pi], [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]))) 
prior = addFactor!(fg, [:x0], PriorPose3( MvNormal([0., 0, 0, 0, 0, pi], [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]))) 

# Drive around in a hexagon vertically
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose3)
    pp = Pose3Pose3(MvNormal([10.0, 0, 0, 0, pi/3, 0], [0.5, 0.5, 0.5, 0.05, 0.05, 0.05]))
    addFactor!(fg, [psym;nsym], pp )
end

ensureAllInitialized!(fg)
getPPESuggested.(fg, sortDFG(ls(fg)))

# tree = solveTree!(fg)


# end
