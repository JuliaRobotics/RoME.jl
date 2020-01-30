
# circle trajectories to show marginalization
using Distributed
addprocs(4)

# using Revise

using RoME, RoMEPlotting
@everywhere using RoME

Gadfly.set_default_plot_size(35cm,25cm)


SIZE = 10

### SHOW marginalization--------------------------------------------------------


# drive first half
fg = generateCanonicalFG_Circle(SIZE, stopEarly=round(Int, SIZE/2), biasTurn=-0.05)

defaultFixedLagOnTree!(fg, round(Int, SIZE/2))

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
# getSolverParams(fg).N = 200
# tree = wipeBuildNewTree!(fg)
# drawTree(tree, filepath=joinLogPath(fg, "bt5.pdf"))


tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg))

plfl1 = drawPosesLandms(fg, spscale=1.0)


generateCanonicalFG_Circle(SIZE, fg=fg, biasTurn=-0.05, loopClosure=true)
# drawGraph(fg, show=true)


ensureAllInitialized!(fg)


plfl2 = drawPosesLandms(fg, spscale=1.0)

# isMarginalized(fg, :x10)
# isSolved(fg, :x10)
# isSolved(fg, :x11)


fg_ = deepcopy(fg)

##

vo = getEliminationOrder(fg)
filter!(x->!(x in [:x5;:x4]), vo)
push!(vo, :x4)
push!(vo, :x5)

tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg), variableOrder=vo)


plfl3 = drawPosesLandms(fg, spscale=1.0)


## draw compound image

plfl4 = drawPosesLandms(fg, spscale=1.0, manualColor="gray30", point_size=4pt)

plfl2



## export plots
plfl1 |> PDF(joinLogPath(fg, "circ$(SIZE)_fl1.pdf") )
plfl2 |> PDF(joinLogPath(fg, "circ$(SIZE)_fl2.pdf") )
plfl3 |> PDF(joinLogPath(fg, "circ$(SIZE)_fl3.pdf") )

saveDFG(fg_, joinLogPath(fg,"fg_before_2nd"))
saveDFG(fg, joinLogPath(fg,"fg_after_2nd"))


### SHOW Incremental Recycling--------------------------------------------------------



# drive first half
fg = generateCanonicalFG_Circle(SIZE, stopEarly=round(Int, SIZE/2)+1, biasTurn=-0.05)


getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true

tree, smt, hist = solveTree!(fg)


plinc1 = drawPosesLandms(fg, spscale=1.0)


## drive second part
generateCanonicalFG_Circle(SIZE, fg=fg, biasTurn=-0.05, loopClosure=true)
# drawGraph(fg, show=true)


ensureAllInitialized!(fg)

plinc2 = drawPosesLandms(fg, spscale=1.0)



fg_ = deepcopy(fg)

tree, smt, hist = solveTree!(fg, tree, recordcliqs=ls(fg))


plinc3 = drawPosesLandms(fg, spscale=1.0)





## export plots
plinc1 |> PDF(joinLogPath(fg, "circ$(SIZE)_inc1.pdf") )
plinc2 |> PDF(joinLogPath(fg, "circ$(SIZE)_inc2.pdf") )
plinc3 |> PDF(joinLogPath(fg, "circ$(SIZE)_inc3.pdf") )

saveDFG(fg_, joinLogPath(fg,"fg_before_2nd"))
saveDFG(fg, joinLogPath(fg,"fg_after_2nd"))




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
