## Caesar.jl, need an easy way to accumulate ground truth from factors only (like rigid transform tree)
# e.g., PriorX0 + Odometry + Odometry

# using Revise

using RoME
using DistributedFactorGraphs

using Test


@testset "Test parametric mean accumulation..." begin

fg = initfg()

addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0;], PriorPose2(MvNormal(zeros(3),0.001*Matrix(LinearAlgebra.I, 3,3))))

addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal([10;0;0.0],0.001*Matrix(LinearAlgebra.I, 3,3))))

# drawGraph(fg)


# add parametric means
val = accumulateFactorMeans(fg, [:x0f1; :x0x1f1])

@test norm(val - [10;0;0]) < 1e-5


end



#
