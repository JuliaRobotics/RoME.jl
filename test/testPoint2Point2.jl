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


@testset "test Point2Point2Range{T}..." begin

N=200
fg = initfg()

addNode!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], eye(2))))

addNode!(fg, :x1, Point2)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], eye(2))))

addNode!(fg, :l1, Point2)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)))
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)))

tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)

@test sum( 90 .< getVal(fg, :l1)[1,:] .< 110 ) > 40
@test sum( -10 .< getVal(fg, :l1)[2,:] .< 10 ) > 40

@test sum( -10 .< getVal(fg, :l1)[1,:] .< 10 ) > 40
@test sum( 90 .< getVal(fg, :l1)[2,:] .< 110 ) > 40

end

# debugging
#
# using RoMEPlotting
#
# plotKDE(fg, [:x0; :x1; :l1], levels=3)
# plotKDE(fg, :l1)
#
# const IIF = IncrementalInference
# const KDE = KernelDensityEstimate
# stuff = IncrementalInference.localProduct(fg, :l1)
# stuff
#
# plotKDE(stuff[2])
#
# pts = IIF.approxConv(fg, :x0l1f1, :l1)
#
# plotKDE(KDE.kde!(pts))

#
