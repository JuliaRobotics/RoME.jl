
# using Revise

using Test
using RoME
using DistributedFactorGraphs
using TensorCast
##

@testset "sanity check on Hex example" begin

##

fg = generateGraph_Hexagonal(graphinit=false)


getSolverParams(fg).drawtree = false
getSolverParams(fg).showtree = false
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).downsolve = true
getSolverParams(fg).multiproc = false
getSolverParams(fg).async = false

# @test_skip false
getSolverParams(fg).useMsgLikelihoods=true
# @error "Restore useMsgLikelihoods=true"

# direct solve would be
tree = solveTree!(fg)
# tree = solveTree!(fg)


@error "Hex init test degraded quality during useMsgLikelihoods refactor, must restore to 55"
M = getManifold(Pose2())

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x0)))[i][j]
@test 35 < sum(-3.0 .< cs[1,:] .< 3.0)
@test 35 < sum(-3.0 .< cs[2,:] .< 3.0)
@test 35 < sum(-0.3 .< cs[3,:] .< 0.3)

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x1)))[i][j]
@test 35 < sum(7.0 .< cs[1,:] .< 13.0)
@test 35 < sum(-3.0 .< cs[2,:] .< 3.0)
@test 35 < sum(0.7 .< cs[3,:] .< 1.3)

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x2)))[i][j]
@test 35 < sum(12.0 .< cs[1,:] .< 18.0)
@test 35 < sum(6.0 .< cs[2,:] .< 11.0)
@test 35 < sum(1.8 .< cs[3,:] .< 2.4)

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x3)))[i][j]
@test 35 < sum(7.0 .< cs[1,:] .< 13.0)
@test 35 < sum(15.0 .< cs[2,:] .< 20.0)
# @test 35 < sum(-0.3 .< cs[3,:] .< 0.3)
μ = mean(M, getPoints(getBelief(fg, :x3)))
@test isapprox(μ.parts[1], [11; 17.5], atol=3.0)
@test isapprox(SpecialOrthogonal(2), μ.parts[2], [-1 0; 0 -1], atol=0.5)
Σ = cov(getManifold(Pose2()), getPoints(getBelief(fg, :x3)))
@test all(diag(Σ) .< [5,5,1].^2) #TODO smaller

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x4)))[i][j]
@test 35 < sum(-5.0 .< cs[1,:] .< 5.0)
@test 35 < sum(13.0 .< cs[2,:] .< 22.0)
@test 35 < sum(-2.8 .< cs[3,:] .< -1.5)

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x5)))[i][j]
@test 35 < sum(-8.0 .< cs[1,:] .< -2.0)
@test 35 < sum(6.0 .< cs[2,:] .< 11.0)
@test 35 < sum(-1.3 .< cs[3,:] .< -0.7)

@cast cs[j,i] := getCoordinates.(Pose2, getPoints(getBelief(fg, :x6)))[i][j]
@test 35 < sum(-3.0 .< cs[1,:] .< 3.0)
@test 35 < sum(-3.0 .< cs[2,:] .< 3.0)
@test 35 < sum(-0.3 .< cs[3,:] .< 0.3)

@cast cs[j,i] := getCoordinates.(Point2, getPoints(getBelief(fg, :l1)))[i][j]
@test 35 < sum(17.0 .< cs[1,:] .< 23.0)
@test 35 < sum(-5.0 .< cs[2,:] .< 5.0)


##

0
end # testset


#
# #  Do some plotting
# using RoMEPlotting
# #
# Gadfly.set_default_plot_size(35cm,20cm)
# pl = drawPosesLandms(fg, drawhist=true)



## OR manually do clique by clique
#
# # solve by hand, one cliq at a time
# tree = wipeBuildNewTree!(fg)
# # drawTree(tree,show=true)
#
#
# # solve the first cliq
# smt, hist = solveCliq!(fg, tree, :x0, recordcliq=true)
#
# # solve second clq
# smt, hist = solveCliq!(fg, tree, :x2, cliqHistories=hist, recordcliq=true)
# # resetCliqSolve!(fg, tree, :x2)
#
# # solve forth clq
# smt, hist = solveCliq!(fg, tree, :x4, cliqHistories=hist, recordcliq=true)
# # solve forth clq
# smt, hist = solveCliq!(fg, tree, :x6, cliqHistories=hist, recordcliq=true)
# # solve forth clq
# smt, hist = solveCliq!(fg, tree, :x3, cliqHistories=hist, recordcliq=true)
