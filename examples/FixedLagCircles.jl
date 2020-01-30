
# circle trajectories to show marginalization
using Distributed
addprocs(4)

# using Revise

using RoME, RoMEPlotting
@everywhere using RoME

Gadfly.set_default_plot_size(35cm,25cm)



### SHOW marginalization--------------------------------------------------------


# drive first half
fg = generateCanonicalFG_Circle(20, stopEarly=10, biasTurn=-0.05)

defaultFixedLagOnTree!(fg, 10)

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
# getSolverParams(fg).N = 200

tree, smt, hist = solveTree!(fg)

drawPosesLandms(fg, spscale=1.0)


generateCanonicalFG_Circle(20, fg=fg, biasTurn=-0.05, loopClosure=true)
# drawGraph(fg, show=true)


ensureAllInitialized!(fg)


drawPosesLandms(fg, spscale=1.0)


isMarginalized(fg, :x10)
isSolved(fg, :x10)
isSolved(fg, :x11)

0




##
fg_ = deepcopy(fg)
tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg))
# fetchAssignTaskHistoryAll!(tree, smt) # if async


pl = drawPosesLandms(fg, spscale=1.0)
pl |> PDF(joinpath(ENV["HOME"], "Downloads", "circ20.pdf"))


### SHOW Recycling--------------------------------------------------------





# drive first half
fg = generateCanonicalFG_Circle(20, stopEarly=10, biasTurn=-0.05)


getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true

tree, smt, hist = solveTree!(fg)



generateCanonicalFG_Circle(20, fg=fg, biasTurn=-0.05, loopClosure=true)
# drawGraph(fg, show=true)


ensureAllInitialized!(fg)

drawPosesLandms(fg, spscale=1.0)



fg_ = deepcopy(fg)
tree, smt, hist = solveTree!(fg, tree, recordcliqs=ls(fg))


drawPosesLandms(fg, spscale=1.0)


## debugging solves


using Gadfly
drawTree(tree, imgs=false, show=true)


spyCliqMat(tree, :x10)


getUpMsgs(tree, :x9)

plotTreeProductUp(fg,tree,:x11)
plotTreeProductUp(fg, tree, :x10, :x11)


# look at plotting with Makie


## new experimental Makie plotting features
using Caesar

include(joinpath((Base.pathof(Caesar) |> splitpath)[1:end-2]..., "examples","marine","auv","Sandshark", "MakiePlotsFG.jl") )


plotVariableBeliefs(fg, r"x\d", resolution=(1920,1080), fade=5, sortVars=true)



0



## ideas for IIF


## check clique messages separately
# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)

# getSolverParams(fg).downsolve = false
# getSolverParams(fg).async = true
# stuffx9 = solveCliq!(fg, tree, :x9, recordcliq=true)
# drawTree(tree)

# umx9 = getUpMsgs(tree,:x9)
# plotKDE(getKDE(fg, :x10), dims=[1;2])
# plotKDE(umx9[:x10][1], dims=[1;2])
#
# stuffx10 = solveCliq!(fg, tree, :x10, recordcliq=true)
#
# umx10 = getUpMsgs(tree,:x10)
# plotKDE(getKDE(fg, :x11), dims=[1;2])
# plotKDE(umx10[:x11][1], dims=[1;2])
#


#
