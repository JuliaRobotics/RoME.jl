using RoME, Distributions

using Base: Test

@testset "basic Point2Point2 test" begin

fg = initfg()

addNode!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), eye(2))))

addNode!(fg, :x1, Point2)
addFactor!(fg, [:x0;:x1], Point2Point2(MvNormal([10;0.0], eye(2))))

tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)

@test sum( abs.(Base.mean(getVal(fg, :x0),2) - [0.0;0]) .< [0.5;0.5]) == 2
@test sum( abs.(Base.mean(getVal(fg, :x1),2) - [10.0;0]) .< [0.5;0.5]) == 2

end
