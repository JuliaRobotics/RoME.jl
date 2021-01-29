using RoME
using Test
using Statistics

## MWE Pose2Point2 from #388
@testset "basic Pose2Point2 test" begin

    fg = initfg()

    addVariable!(fg, :x1, Pose2)
    addVariable!(fg, :l1, Point2)

    addFactor!(fg, [:x1], PriorPose2(MvNormal([0.,0, 0], [0.01, 0.01, 0.01])))

    addFactor!(fg, [:x1; :l1], Pose2Point2(MvNormal([0.0,-1], [0.1,0.1])))

    ensureAllInitialized!(fg)

    tree, smt, hist = solveTree!(fg)

    @test isapprox(mean(getVal(fg, :x1),dims=2), [0,0,0], atol = 0.05)
    @test isapprox(mean(getVal(fg, :l1),dims=2), [0,-1], atol = 0.05)
end


