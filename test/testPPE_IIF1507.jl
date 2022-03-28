
using RoME
using Test


##

@testset "Wrong factor type behavior, IIF 1507" begin
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

# purposefully send "wrong" factor type
# see discussion: https://github.com/JuliaRobotics/IncrementalInference.jl/issues/1507#issuecomment-1081087381
addFactor!(fg, [:x0;:l0], Point2Point2Range(Normal(5, 1)) )


##



##
end


##
