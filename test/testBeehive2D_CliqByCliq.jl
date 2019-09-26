
# using Revise

using Test

using RoME


include(joinpath(@__DIR__,"BeehiveTestUtils.jl"))


@testset "sanity check on Hex example" begin

## start with an empty factor graph object
fg = initfg()

posecount = 0

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)
posecount += 1


# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([0.0; 0.0; 0.0],
                                           Matrix(Diagonal([0.1;0.1;0.05].^2))) ), autoinit=false )

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, labels=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false )



## hex 1

posecount = driveHex(fg, posecount)

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg,engine="neato")

getSolverParams(fg).drawtree = false
getSolverParams(fg).showtree = false
getSolverParams(fg).downsolve = false
getSolverParams(fg).multiproc = false
getSolverParams(fg).async = false
# to disable parent factor sharing
# getSolverParams(fg).devParams[:dontUseParentFactorsInitDown] = ""



# direct solve would be
tree, smt, hist = solveTree!(fg) #, recordcliqs=ls(fg))


@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[1,:] .< 3.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[2,:] .< 3.0)
@test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x0))[3,:] .< 0.3)

@test 80 < sum(7.0 .< getPoints(getKDE(fg, :x1))[1,:] .< 13.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x1))[2,:] .< 3.0)
@test 80 < sum(0.7 .< getPoints(getKDE(fg, :x1))[3,:] .< 1.3)

@test 80 < sum(12.0 .< getPoints(getKDE(fg, :x2))[1,:] .< 18.0)
@test 80 < sum(6.0 .< getPoints(getKDE(fg, :x2))[2,:] .< 11.0)
@test 80 < sum(1.8 .< getPoints(getKDE(fg, :x2))[3,:] .< 2.4)

@test 80 < sum(7.0 .< getPoints(getKDE(fg, :x3))[1,:] .< 13.0)
@test 80 < sum(15.0 .< getPoints(getKDE(fg, :x3))[2,:] .< 20.0)
# @test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x3))[3,:] .< 0.3)

@test 80 < sum(-5.0 .< getPoints(getKDE(fg, :x4))[1,:] .< 5.0)
@test 80 < sum(13.0 .< getPoints(getKDE(fg, :x4))[2,:] .< 22.0)
@test 80 < sum(-2.8 .< getPoints(getKDE(fg, :x4))[3,:] .< -1.5)

@test 80 < sum(-8.0 .< getPoints(getKDE(fg, :x5))[1,:] .< -2.0)
@test 80 < sum(6.0 .< getPoints(getKDE(fg, :x5))[2,:] .< 11.0)
@test 80 < sum(-1.3 .< getPoints(getKDE(fg, :x5))[3,:] .< -0.7)

@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[1,:] .< 3.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[2,:] .< 3.0)
@test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x6))[3,:] .< 0.3)

@test 80 < sum(17.0 .< getPoints(getKDE(fg, :l1))[1,:] .< 23.0)
@test 80 < sum(-5.0 .< getPoints(getKDE(fg, :l1))[2,:] .< 5.0)



# #  Do some plotting
# using RoMEPlotting
#
# Gadfly.set_default_plot_size(35cm,20cm)
# drawPosesLandms(fg, meanmax=:max) |> PDF("/tmp/caesar/test.pdf");  @async run(`evince /tmp/caesar/test.pdf`)
# 0
# # drawTree(tree, imgs=true)


end # testset



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









# ## RANDOM DEBUG DEV CODE BELOW
#
# hist = getCliqSolveHistory(tree, :x1)
#
# tree, smt, hist = solveTree!(fg, skipcliqids=[:x1;:x6;:x4;:x3], recordcliqs=[:x2;:x0])
#
# plotTreeUpMsgs(fg, tree, :x1, levels=1)
#
#
# cliq = getCliq(tree, :x2)
#
# # OLD
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
# plotKDE(dwinmsgs[:x1][1], dims=[1;2], levels=2)
#
# # NEW
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
# plotKDE(dwinmsgs[:x1][1], dims=[1;2], levels=2)
#
#
#
# ## Check init message for x3
#
# cliq = getCliq(tree, :x1)
#
# # OLD
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
# # NEW
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
#
#
# ## Check init message for x4
#
# cliq = getCliq(tree, :x4)
#
# # OLD
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
# # NEW
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
#
#
#
# plotKDE(fg, :x1, dims=[1;2])
# plotPose(fg, :x1)
#
# plotKDE(fg, [:x0;:x1;:x2;:x4;:x6], dims=[1;2],levels=1)
#
#
# drawTree(tree, imgs=true)
#
#
#
# ## check contents
#
# cliq = getCliq(tree, :x1)
#
# getCliqMsgsUp(cliq)
#
#
# #
#
# # DEBUG after solveCliq :x2
#
#
# # solve third clq
# smt, hist = solveCliq!(fg, tree, :x1, cliqHistories=hist, recordcliq=true)
#
#
# # Plot to see what is going on
# Gadfly.set_default_plot_size(35cm,20cm)
# plotKDE(fg, [:x0;:x1;:x2], dims=[1;2],levels=1)
#
# plotKDE(fg, :x1, dims=[1;2],levels=1)
#
#
#
# ## DEBUG where is down message :x3
#
# prnt = getCliq(tree, :x3)
# getCliqInitUpMsgs(prnt)
#
#
# assignTreeHistory!(tree, hist)
# printCliqHistorySummary(tree, :x1)
#
#
# # think the issue is in here
# # 21:01:00.856  8   null         attemptCliqInitUp     false null | upsolved upsolved
#
#
#
# # somethings up with cliq x1
#
# csmc_8_test = getCliqSolveHistory(tree, :x1)[8][4]
# csmc_9_test = getCliqSolveHistory(tree, :x1)[9][4]
#
# # Just before the mistake
# drawTree(csmc_8_test.tree, show=true)
# prnt_8_test = getCliq(csmc_8_test.tree, :x3)
# getCliqInitUpMsgs(prnt_8_test)
#
# # just after the mistake
# prnt_9_test = getCliq(csmc_9_test.tree, :x3)
# getCliqInitUpMsgs(prnt_9_test)
#
#
# # so develop  based on sandbox step 8
# stuff = sandboxCliqResolveStep(tree,:x1,8)
# getCliqInitUpMsgs(getCliq(stuff[4].tree,:x3))
#
#
#
# IIF.getCliqMsgsUp(tree,:x1)
#
#
# ## DEBUG
#
#
#
#
#
#
#
# ## DEBUG after solveCliq :x2 difference in down init cycle when :dontUseParentFactorsInitDown
#
# assignTreeHistory!(tree, hist)
# printCliqHistorySummary(tree, :x2)
#
# # WORKS
#
# delete!(getSolverParams(getCliqSolveHistory(tree, :x2)[15][4].dfg).devParams, :dontUseParentFactorsInitDown)
#
# good = sandboxCliqResolveStep(tree,:x2,15)
#
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x1, title="(X2,) sandbox 15 good,")
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x2, title="(X2,) sandbox 15 good,")
#
#
# # BREAKS
#
# getSolverParams(getCliqSolveHistory(tree, :x2)[15][4].dfg).devParams[:dontUseParentFactorsInitDown] = ""
#
# bad = sandboxCliqResolveStep(tree,:x2,15)
#
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x1, title="(X2,) sandbox 15 bad,")
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x2, title="(X2,) sandbox 15 bad,")
#
# plotPairVariables(good[4].cliqSubFg, bad[4].cliqSubFg, :x2, title="(X2,) sandbox 15,")
# # plotPairPose2(good[4].cliqSubFg, bad[4].cliqSubFg, :x2, title="(X2,) sandbox 15,")
#
# Gadfly.set_default_plot_size(30cm,18cm)
#
# ## DEBUG
#
#
#
#
#
#
#
# #3 First difference, up msg from (:x2) cliq is 'wrong'
#
# Gadfly.set_default_plot_size(30cm,18cm)
# um1 = IIF.getCliqMsgsUp(tree, :x0)
# plotKDE(um1[:x1], dims=[1;2], levels=2)
#
#
#
# um2 = IIF.getCliqMsgsUp(tree, :x2)
# plotKDE(um2[:x1], dims=[1;2], levels=2)
#
#
#
#
# csmcX2_test = getCliqSolveHistory(tree, :x2)
#
# # csmcX2_test[15][4].cliqSubFg
#
# plotKDE(getKDE(csmcX2_test[17][4].cliqSubFg, :x1),dims=[1;2],levels=2,title="csmcX2_test[17]")
#
#
#
#
# ## SOMETHING changes between csmcX2_ref[16/17] and csmcX2_test[16/17]
#
# writeGraphPdf(csmcX2_ref[17][4].cliqSubFg, show=true)
# writeGraphPdf(csmcX2_test[17][4].cliqSubFg, show=true)
#
#
# plotKDE(csmcX2_ref[16][4].cliqSubFg, [:x1;:x2;:x3], dims=[1;2])
# plotKDE(csmcX2_test[16][4].cliqSubFg, [:x1;:x2;:x3], dims=[1;2])
#
#
#
# isInitialized(csmcX2_test[16][4].cliqSubFg, :x3)
#
#
#
#
#
#
# # before init cycle
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x1)
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x2)
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x3)
#
#
# # after solve
# plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x1)
# plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x2)
# plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x3)
