
using RoME
using Test


##

@testset "Test PPE not being populated by solve, IIF 1507" begin
##

fg = initfg()

addVariable!(fg, :x0, Pose2)
addVariable!(fg, :x1, Pose2)
addVariable!(fg, :x2, Pose2)
addVariable!(fg, :l0, Point2)

addFactor!(fg, [:x0], PriorPose2( MvNormal(zeros(3), diagm([0.1, 0.1, 0.1]))  ) )
addFactor!(fg, [:x0; :x1], Pose2Pose2( MvNormal([1, 1, pi / 3], diagm([0.1, 0.1, 0.1])) ) )
addFactor!(fg, [:x1; :x2], Pose2Pose2( MvNormal([1, 1, pi / 3], diagm([0.1, 0.1, 0.1])) ) )
addFactor!(fg, [:l0], PriorPoint2( MvNormal([5, 0], diagm([2, 2])) ) )
addFactor!(fg, [:x0;:l0], Point2Point2Range(Normal(5, 1)) )


##



##
end


##
