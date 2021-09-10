using RoME
using Test
using Manifolds: ProductRepr

@testset "Basic PartialPriorYawPose2 test" begin

fg = initfg()

# fg.solverParams.graphinit=false
addVariable!(fg, :x, Pose2)
addVariable!(fg, :l, Point2)

addFactor!(fg, [:x], PartialPriorYawPose2(MvNormal([0.,0, 0], [0.01, 0.01, 0.01])))
addFactor!(fg, [:l], PriorPoint2(MvNormal([1.,1], [0.01, 0.01])))

addFactor!(fg, [:x; :l], Pose2Point2(MvNormal([1.,1], [0.1,0.1])))

initAll!(fg)

tree = solveTree!(fg)

M = getManifold(Pose2)
@test isapprox(M, mean(M, getVal(fg, :x)), ProductRepr([0,0], [1 0; 0 1]), atol=0.05) 
@test isapprox(mean(getVal(fg, :l)), [1,1], atol = 0.05)
end
