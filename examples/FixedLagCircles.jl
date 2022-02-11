
# circle trajectories to show marginalization
using Distributed
addprocs(4)

using RoME, RoMEPlotting, DistributedFactorGraphs, KernelDensityEstimatePlotting
@everywhere using RoME

Gadfly.set_default_plot_size(35cm,25cm)

include(joinpath(@__DIR__, "FixedLagCirclePlotting.jl"))


SIZE = 10

### SETUP first half of circle----------------------------------------------

# drive first half
fg = generateGraph_Circle(SIZE, stopEarly=round(Int, SIZE/2), biasTurn=-0.05)

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
# getSolverParams(fg).N = 200
# tree = wipeBuildNewTree!(fg)
# drawTree(tree, filepath=joinLogPath(fg, "bt5.pdf"))


tree = solveTree!(fg, recordcliqs=ls(fg))

# plfl1 = drawPosesLandms(fg, spscale=1.0)
fg5a = deepcopy(fg)
tree5 = deepcopy(tree)

plotCirc10BA(fg5a, fg5a, filepath=joinLogPath(fg, "circ5A.pdf"), pointsList_=[:l1])
pl5blue = plotCirc10BA(fg5a, fg5a, filepath=joinLogPath(fg, "circ5A_blue.pdf"), lineColor="deepskyblue")

saveDFG(fg5a, joinLogPath(fg5a,"fg_5A"))

# drive second half

generateGraph_Circle(SIZE, fg=fg, biasTurn=-0.05, loopClosure=true, kappaOdo=3)

initAll!(fg)

# keep copy for later
fg_ = deepcopy(fg)

saveDFG(fg_, joinLogPath(fg_,"fg_10B"))



### SHOW Incremental Recycling--------------------------------------------------------

fgI = deepcopy(fg_)

tree = solveTree!(fgI, deepcopy(tree5), recordcliqs=ls(fgI))


plotCirc10BA(fg_, fgI, filepath=joinLogPath(fg, "circ$(SIZE)IncrBA.pdf"), pointsList_=[:x8])


saveDFG(fgI, joinLogPath(fgI,"fg_IncrA"))


calcCliquesRecycled(tree5)
calcCliquesRecycled(tree)


## SHOW MARGINALIZATION---------------------------------------------------------

fgM = deepcopy(fg_)

defaultFixedLagOnTree!(fgM, round(Int, SIZE/2))

# use slightly special variable ordering
vo = getEliminationOrder(fgM)
filter!(x->!(x in [:x5;:x4]), vo)
push!(vo, :x4)
push!(vo, :x5)

tree = solveTree!(fgM, recordcliqs=ls(fgM), variableOrder=vo)

## Plot ellipses to illustrate covariance fit

pl10marg = plotCirc10BA(fg_, fgM, filepath=joinLogPath(fgM, "circ$(SIZE)MargBA.pdf"), pointsList=[:x5])

for ll in pl10marg.layers
  push!(pl5blue.layers, ll)
end
pl5blue |> PDF(joinLogPath(fgM, "circ$(SIZE)MargBA_blue.pdf"), 10cm, 8cm)



saveDFG(fgM, joinLogPath(fgM,"fg_MargA"))


# need a fix
calcCliquesRecycled(tree)


### DEBUGGING

cliq = getClique(tree, :x1)
data = getData(cliq).allmarginalized


## SHOW BOTH MARG AND INCR TOGETHER---------------------------------------------

fgR = deepcopy(fg_)

defaultFixedLagOnTree!(fgR, round(Int, SIZE/2), limitfixeddown=false)


tree = solveTree!(fgR, deepcopy(tree5), recordcliqs=ls(fgR))


pl10recl = plotCirc10BA(fg_, fgR, filepath=joinLogPath(fg, "circ$(SIZE)ReclBA.pdf"), levels=3, drawEllipse=true, ellipseList_=[:x9])
#, contourList_=[:x8], width=11cm, ellipseColor="gray20"


pl5recl = drawPosesLandms(fgR, spscale=1.5, manualColor="cyan4", point_size=4pt, drawhist=false, contour=false, levels=2, lbls=false, to=5, line_width=2pt)

for ll in pl10recl.layers
  push!(pl5recl.layers, ll)
end
pl5recl |> PDF(joinLogPath(fg, "circ10Recl_cyan4.pdf"), 10cm, 8cm)



saveDFG(fgR, joinLogPath(fgR,"fg_ReclA"))









## draw arcs

plotCirc10BA(fg5a, fg5a, filepath=joinLogPath(fg, "circ5A.pdf"))



tfg = initfg()

# add a starting point (skipping prior for brevity)
addVariable!(tfg, :a, Pose2)
initManual!(tfg, :a, 0.01*randn(3,100))

addVariable!(tfg, :a_drt, Pose2)

drt = MutablePose2Pose2Gaussian(MvNormal(zeros(3),Matrix{Float64}(LinearAlgebra.I,3,3)))
addFactor!(tfg, [:a; :a_drt], drt)

Qc = 0.001*Matrix{Float64}(LinearAlgebra.I,3,3)

N = 10
XYT = zeros(N,3)

for i in 1:N
  accumulateDiscreteLocalFrame!(drt,[0.1;0;0.05],Qc)
  drval = accumulateFactorMeans(tfg, [:aa_drtf1;])
  XYT[i,:] = drval
end

XYT




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
