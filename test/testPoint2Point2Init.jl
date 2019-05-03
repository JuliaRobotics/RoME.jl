
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
# eo = getEliminationOrder(fg, ordering=:qr)
# eo = [1;3;5]
# tree = buildTreeFromOrdering!(fgl,eo)


cliq = tree.cliques[2]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=false, limititers=1 )


cliq = tree.cliques[1]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=false, limititers=1 )


@test 20 < sum( 90 .< getVal(fg, :l1)[1,:] .< 110 )
@test 20 < sum( -10 .< getVal(fg, :l1)[2,:] .< 10 )

@test 20 < sum( -10 .< getVal(fg, :l1)[1,:] .< 10 )
@test 20 < sum( 90 .< getVal(fg, :l1)[2,:] .< 110 )



@test 80 < sum( 90 .< getVal(fg, :x0)[1,:] .< 110 )
@test 80 < sum( -10 .< getVal(fg, :x0)[2,:] .< 10 )


@test 80 < sum( -10 .< getVal(fg, :x1)[1,:] .< 10 )
@test 80 < sum( 90 .< getVal(fg, :x1)[2,:] .< 110 )



ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=true);



@test 20 < sum( 90 .< getVal(fg, :l1)[1,:] .< 110 )
@test 20 < sum( -10 .< getVal(fg, :l1)[2,:] .< 10 )

@test 20 < sum( -10 .< getVal(fg, :l1)[1,:] .< 10 )
@test 20 < sum( 90 .< getVal(fg, :l1)[2,:] .< 110 )



@test 80 < sum( 90 .< getVal(fg, :x0)[1,:] .< 110 )
@test 80 < sum( -10 .< getVal(fg, :x0)[2,:] .< 10 )


@test 80 < sum( -10 .< getVal(fg, :x1)[1,:] .< 10 )
@test 80 < sum( 90 .< getVal(fg, :x1)[2,:] .< 110 )


end
