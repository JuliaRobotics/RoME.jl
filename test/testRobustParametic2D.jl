using RoME
using Test

@testset "Triangulation test in 2D, 4 beacons with outlier" begin
# landmark prior noise 
lmp_noise = diagm([0.01;0.01].^2)

# create a new factor graph and set to not auto initialize
fg = initfg()
fg.solverParams.graphinit = false

# add 3 landmarks variables for inliers points (TranslationGroup(2)) with prior mesurements for each
addVariable!(fg, :l1, Point2)
addVariable!(fg, :l2, Point2)
addVariable!(fg, :l3, Point2)

addFactor!(fg, [:l1], PriorPoint2(MvNormal([10.0;1.0],lmp_noise)))  
addFactor!(fg, [:l2], PriorPoint2(MvNormal([10.0+sqrt(3)/2;-0.5],lmp_noise)))
addFactor!(fg, [:l3], PriorPoint2(MvNormal([10.0-sqrt(3)/2;-0.5],lmp_noise)))

# Add a 2D Pose (SE2) with measurements to landmarks for pose to end up at (10, 0, 0)
addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x1, :l1], Pose2Point2Bearing(Normal(pi/2,0.05)))
addFactor!(fg, [:x1, :l2], Pose2Point2Bearing(Normal(-pi/6,0.05)))
addFactor!(fg, [:x1, :l3], Pose2Point2Bearing(Normal(-pi+pi/6,0.05)))

#solve the graph - calls solve_RLM which uses Manopt.LevenbergMarquardt!
IIF.solveGraphParametric!(fg)
# get the result of the pose
x1_in = getVal(fg, :x1; solveKey=:parametric)[1]

@test isapprox(x1_in, ArrayPartition([10.0, 0], [1.0 0; 0 1.0]))

# Add the outlier landmark (l4) and measurements
addVariable!(fg, :l4, Point2)
addFactor!(fg, [:l4], PriorPoint2(MvNormal([10.5;1.0],lmp_noise)))  
addFactor!(fg, [:x1, :l4], Pose2Point2Bearing(Normal(pi/2,0.05)))

IIF.solveGraphParametric!(fg)
x1_out = getVal(fg, :x1; solveKey=:parametric)[1]
# the outlier scews the result to 
# ([10.163, -0.068], [0.9970 0.077; -0.077 0.997])

end