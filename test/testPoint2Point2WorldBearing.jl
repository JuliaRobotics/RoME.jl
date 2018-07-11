# testing Point2Point2WorldBearing factor

using RoME, IncrementalInference, Distributions

using Base: Test




@testset "test Point2Point2WorldBearing factor is working" begin

N = 100
fg = initfg()

addNode!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), eye(2))))

ensureAllInitialized!(fg)

addNode!(fg, :x1, Point2)
# addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;-30.0], 100*eye(2))))


pp2 = Point2Point2WorldBearing(Normal(-3pi/4,0.1))
addFactor!(fg, [:x0;:x1], pp2)


addNode!(fg, :x2, Point2)
addFactor!(fg, [:x2], PriorPoint2(MvNormal([0.0;-100.0], eye(2))))


pp2 = Point2Point2WorldBearing(Normal(3pi/4,0.1))
addFactor!(fg, [:x2;:x1], pp2)


tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)


pts = getVal(fg, :x1)

@test 0.8*N < sum(-75 .< pts[1,:] .< -25)
@test 0.8*N < sum(-75 .< pts[2,:] .< -25)


@test 0 <= sum(0 .< pts[1,:]) < 0.05*N
@test 0 <= sum(0 .< pts[2,:]) < 0.05*N


@test 0 <= sum(pts[1,:] .< -100.0) < 0.05*N
@test 0 <= sum(pts[2,:] .< -100.0) < 0.05*N


end


# writeGraphPdf(fg)

#
# using RoMEPlotting, Gadfly
# xmin=-150
# xmax=50
# ymin=-150
# ymax=50
#
#
#
# pl = plotKDE(fg, [:x0, :x1], dims=[1;2]);
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl
#
# pl = plotKDE(fg, [:x0, :x1, :x2], dims=[1;2]);
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl





#
# gg = (x) -> pdf(Beta(2.0,6.0), x)
#
# g2 = (x) -> pdf(Rayleigh(100), x)
#
#
# gg(0.5)
#
# plot(g2, -5.0,250.0)
#



#
