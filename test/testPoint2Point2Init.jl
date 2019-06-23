
using RoME

using Test



@testset "test tree based autoinit on small bimodal point2 trilateration..." begin


N=100
fg = initfg()

addVariable!(fg, :x0, Point2, N=N)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], Matrix{Float64}(LinearAlgebra.I, 2,2))), autoinit=false)

addVariable!(fg, :x1, Point2, N=N)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], Matrix{Float64}(LinearAlgebra.I, 2,2))), autoinit=false)

addVariable!(fg, :l1, Point2, N=N)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)), autoinit=false)
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)), autoinit=false)


tree = wipeBuildNewTree!(fg)
# drawTree(tree, filepath="/tmp/caesar/bt.pdf", show=true)

# eo = getEliminationOrder(fg, ordering=:qr)
# eo = [1;3;5]
# tree = buildTreeFromOrdering!(fg,eo)

batchSolve!(fg)

# cliq = tree.cliques[2]
# clst = cliqInitSolveUpByStateMachine!(fg, tree, cliq, drawtree=false, limititers=15 )

# cliq = tree.cliques[1]
# clst = cliqInitSolveUpByStateMachine!(fg, tree, cliq, drawtree=false, limititers=15 )


# using RoMEPlotting
#
# plotKDE(fg, [:x0;:x1]) |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)
# plotKDE(fg, :l1) |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)


@test 20 < sum( 80 .< getVal(fg, :l1)[1,:] .< 120 )
@test 20 < sum( -20 .< getVal(fg, :l1)[2,:] .< 20 )

@test 20 < sum( -20 .< getVal(fg, :l1)[1,:] .< 20 )
@test 20 < sum( 80 .< getVal(fg, :l1)[2,:] .< 120 )



@test 80 < sum( 80 .< getVal(fg, :x0)[1,:] .< 120 )
@test 80 < sum( -20 .< getVal(fg, :x0)[2,:] .< 20 )


@test 80 < sum( -20 .< getVal(fg, :x1)[1,:] .< 20 )
@test 80 < sum( 80 .< getVal(fg, :x1)[2,:] .< 120 )



# ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
# downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=false);



@test 20 < sum( 80 .< getVal(fg, :l1)[1,:] .< 120 )
@test 20 < sum( -20 .< getVal(fg, :l1)[2,:] .< 20 )

@test 20 < sum( -20 .< getVal(fg, :l1)[1,:] .< 20 )
@test 20 < sum( 80 .< getVal(fg, :l1)[2,:] .< 120 )



@test 80 < sum( 80 .< getVal(fg, :x0)[1,:] .< 120 )
@test 80 < sum( -20 .< getVal(fg, :x0)[2,:] .< 20 )


@test 80 < sum( -20 .< getVal(fg, :x1)[1,:] .< 20 )
@test 80 < sum( 80 .< getVal(fg, :x1)[2,:] .< 120 )


end
