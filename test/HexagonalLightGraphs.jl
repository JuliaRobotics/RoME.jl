using RoME
using DistributedFactorGraphs
using Test

@testset "Test Hexagonal specifically with LightGraphsDFG..." begin

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



@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[1,:] .< 3.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[2,:] .< 3.0)
@test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x0))[3,:] .< 0.3)

@test 80 < sum(7.0 .< getPoints(getKDE(fg, :x1))[1,:] .< 13.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x1))[2,:] .< 3.0)
@test 80 < sum(0.7 .< getPoints(getKDE(fg, :x1))[3,:] .< 1.3)

@test 80 < sum(12.0 .< getPoints(getKDE(fg, :x2))[1,:] .< 18.0)
@test 80 < sum(6.0 .< getPoints(getKDE(fg, :x2))[2,:] .< 11.0)
@test 80 < sum(1.8 .< getPoints(getKDE(fg, :x2))[3,:] .< 2.4)

@test 80 < sum(7.0 .< getPoints(getKDE(fg, :x3))[1,:] .< 13.0)
@test 80 < sum(15.0 .< getPoints(getKDE(fg, :x3))[2,:] .< 20.0)
# @test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x3))[3,:] .< 0.3)

@test 80 < sum(-4.0 .< getPoints(getKDE(fg, :x4))[1,:] .< 4.0)
@test 80 < sum(15.0 .< getPoints(getKDE(fg, :x4))[2,:] .< 20.0)
@test 80 < sum(-2.4 .< getPoints(getKDE(fg, :x4))[3,:] .< -1.8)

@test 80 < sum(-8.0 .< getPoints(getKDE(fg, :x5))[1,:] .< -2.0)
@test 80 < sum(6.0 .< getPoints(getKDE(fg, :x5))[2,:] .< 11.0)
@test 80 < sum(-1.3 .< getPoints(getKDE(fg, :x5))[3,:] .< -0.7)

@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[1,:] .< 3.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[2,:] .< 3.0)
@test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x6))[3,:] .< 0.3)

@test 80 < sum(17.0 .< getPoints(getKDE(fg, :l1))[1,:] .< 23.0)
@test 80 < sum(-5.0 .< getPoints(getKDE(fg, :l1))[2,:] .< 5.0)


end
