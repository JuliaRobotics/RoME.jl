
using RoME

using Test

##

@testset "test tree based autoinit on small bimodal point2 trilateration..." begin

##

N=100
fg = initfg()
getSolverParams(fg).inflation=250.0

addVariable!(fg, :x0, Point2, N=N)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], Matrix{Float64}(LinearAlgebra.I, 2,2))), graphinit=false)

addVariable!(fg, :x1, Point2, N=N)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], Matrix{Float64}(LinearAlgebra.I, 2,2))), graphinit=false)

addVariable!(fg, :l1, Point2, N=N)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)), graphinit=false)
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)), graphinit=false)

##

tree, smt, hist = solveTree!(fg)

##

# using RoMEPlotting
# Gadfly.set_default_plot_size(25cm,20cm)

# plotKDE(fg, [:x0;:x1]) # |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)
# plotKDE(fg, :l1) # |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)

##

IIF._getCCW(fg, :x0l1f1).inflation = 1000.0
pts = approxConv(fg, :x0l1f1, :l1)

# plotKDE(manikde!(pts, Point2))


##

@test 10 < sum( 80 .< getVal(fg, :l1)[1,:] .< 120 )
@test 10 < sum( -20 .< getVal(fg, :l1)[2,:] .< 20 )

@test 10 < sum( -20 .< getVal(fg, :l1)[1,:] .< 20 )
@test 10 < sum( 80 .< getVal(fg, :l1)[2,:] .< 120 )



@test 80 < sum( 80 .< getVal(fg, :x0)[1,:] .< 120 )
@test 80 < sum( -20 .< getVal(fg, :x0)[2,:] .< 20 )


@test 80 < sum( -20 .< getVal(fg, :x1)[1,:] .< 20 )
@test 80 < sum( 80 .< getVal(fg, :x1)[2,:] .< 120 )




@test 10 < sum( 80 .< getVal(fg, :l1)[1,:] .< 120 )
@test 10 < sum( -20 .< getVal(fg, :l1)[2,:] .< 20 )

@test 10 < sum( -20 .< getVal(fg, :l1)[1,:] .< 20 )
@test 10 < sum( 80 .< getVal(fg, :l1)[2,:] .< 120 )



@test 80 < sum( 80 .< getVal(fg, :x0)[1,:] .< 120 )
@test 80 < sum( -20 .< getVal(fg, :x0)[2,:] .< 20 )


@test 80 < sum( -20 .< getVal(fg, :x1)[1,:] .< 20 )
@test 80 < sum( 80 .< getVal(fg, :x1)[2,:] .< 120 )


end

##