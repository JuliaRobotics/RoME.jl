# testing Point2Point2WorldBearing factor

using RoME
# , IncrementalInference, Distributions

using Test




@testset "test Point2Point2WorldBearing factor is working" begin

global N = 100
global fg = initfg()

addVariable!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), Matrix{Float64}(LinearAlgebra.I, 2,2))))

ensureAllInitialized!(fg)

addVariable!(fg, :x1, Point2)
# addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;-30.0], 100*Matrix{Float64}(LinearAlgebra.I, 2,2))))


global pp2 = Point2Point2WorldBearing(Normal(-3pi/4,0.1))
addFactor!(fg, [:x0;:x1], pp2)


addVariable!(fg, :x2, Point2)
addFactor!(fg, [:x2], PriorPoint2(MvNormal([0.0;-100.0], Matrix{Float64}(LinearAlgebra.I, 2,2))))


global pp2 = Point2Point2WorldBearing(Normal(3pi/4,0.1))
addFactor!(fg, [:x2;:x1], pp2)


tree, smt, hist = solveTree!(fg)


global pts = getVal(fg, :x1)

@test 0.5*N < sum(-70 .< pts[1,:] .< -30)
@test 0.5*N < sum(-70 .< pts[2,:] .< -30)


@test 0 <= sum(0 .< pts[1,:]) < 0.1*N
@test 0 <= sum(0 .< pts[2,:]) < 0.1*N


@test 0 <= sum(pts[1,:] .< -100.0) < 0.15*N
@test 0 <= sum(pts[2,:] .< -100.0) < 0.15*N


end

#
# drawGraph(fg, show = true)
# drawTree(tree, show = true)
#
#
#
#
# #
# using RoMEPlotting, Gadfly
# Gadfly.set_default_plot_size(35cm, 20cm)
# xmin=-150
# xmax=50
# ymin=-150
# ymax=50
#
#
#
# pl = plotKDE(fg, [:x1], dims=[1;2]);
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl
#
# pl = plotKDE(fg, [:x0, :x1], dims=[1;2]);
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl
#
# pl = plotKDE(fg, [:x0, :x1, :x2], dims=[1;2]);
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl
# #
# #
# # drawGraph(fg)
#
# pl = plotLocalProduct(fg, :x1)
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl
#
#
# pl = plotKDE(upMsg(tree,:x0)[:x1][1], levels=1)
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl
#
#
# pts = treeProductUp(fg, tree, :x1, :x1)[1]
# pl = plotKDE(manikde!(pts, Point2), levels=1)
# pl.coord = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# pl

#
# 0

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
