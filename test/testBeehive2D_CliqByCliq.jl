
# using Revise

using Test

using RoME



@testset "sanity check on Hex example" begin

# graphinit=false, regression since DFG v0.6.0 IIF v0.9.0 -- on tree init, issue IIF#602
fg = generateCanonicalFG_Hexagonal(graphinit=false)


getSolverParams(fg).drawtree = false
getSolverParams(fg).showtree = false
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).downsolve = true
getSolverParams(fg).multiproc = false
getSolverParams(fg).async = false
getSolverParams(fg).useMsgLikelihoods=true
# to disable parent factor sharing
# getSolverParams(fg).devParams[:dontUseParentFactorsInitDown] = ""


# direct solve would be
# getSolverParams(fg).async = true ## MANUAL TEMP ONLY
tree, smt, hist = solveTree!(fg) #, recordcliqs=ls(fg))
tree, smt, hist = solveTree!(fg) #, recordcliqs=ls(fg))

# tree.cliques[3]
# drawTree(tree, show=true)
# smt[3]

@test 55 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[1,:] .< 3.0)
@test 55 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[2,:] .< 3.0)
@test 55 < sum(-0.3 .< getPoints(getKDE(fg, :x0))[3,:] .< 0.3)

@test 55 < sum(7.0 .< getPoints(getKDE(fg, :x1))[1,:] .< 13.0)
@test 55 < sum(-3.0 .< getPoints(getKDE(fg, :x1))[2,:] .< 3.0)
@test 55 < sum(0.7 .< getPoints(getKDE(fg, :x1))[3,:] .< 1.3)

@test 55 < sum(12.0 .< getPoints(getKDE(fg, :x2))[1,:] .< 18.0)
@test 55 < sum(6.0 .< getPoints(getKDE(fg, :x2))[2,:] .< 11.0)
@test 55 < sum(1.8 .< getPoints(getKDE(fg, :x2))[3,:] .< 2.4)

@test 55 < sum(7.0 .< getPoints(getKDE(fg, :x3))[1,:] .< 13.0)
@test 55 < sum(15.0 .< getPoints(getKDE(fg, :x3))[2,:] .< 20.0)
# @test 55 < sum(-0.3 .< getPoints(getKDE(fg, :x3))[3,:] .< 0.3)

@test 55 < sum(-5.0 .< getPoints(getKDE(fg, :x4))[1,:] .< 5.0)
@test 55 < sum(13.0 .< getPoints(getKDE(fg, :x4))[2,:] .< 22.0)
@test 55 < sum(-2.8 .< getPoints(getKDE(fg, :x4))[3,:] .< -1.5)

@test 55 < sum(-8.0 .< getPoints(getKDE(fg, :x5))[1,:] .< -2.0)
@test 55 < sum(6.0 .< getPoints(getKDE(fg, :x5))[2,:] .< 11.0)
@test 55 < sum(-1.3 .< getPoints(getKDE(fg, :x5))[3,:] .< -0.7)

@test 55 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[1,:] .< 3.0)
@test 55 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[2,:] .< 3.0)
@test 55 < sum(-0.3 .< getPoints(getKDE(fg, :x6))[3,:] .< 0.3)

@test 55 < sum(17.0 .< getPoints(getKDE(fg, :l1))[1,:] .< 23.0)
@test 55 < sum(-5.0 .< getPoints(getKDE(fg, :l1))[2,:] .< 5.0)




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
