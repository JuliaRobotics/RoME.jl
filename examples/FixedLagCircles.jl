
# circle trajectories to show marginalization
using Distributed
addprocs(4)

using RoME, RoMEPlotting, DistributedFactorGraphs
@everywhere using RoME

Gadfly.set_default_plot_size(35cm,25cm)

include(joinpath(@__DIR__, "FixedLagCirclePlotting.jl"))


SIZE = 10

### SETUP first half of circle----------------------------------------------

# drive first half
fg = generateCanonicalFG_Circle(SIZE, stopEarly=round(Int, SIZE/2), biasTurn=-0.05)

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
# getSolverParams(fg).N = 200
# tree = wipeBuildNewTree!(fg)
# drawTree(tree, filepath=joinLogPath(fg, "bt5.pdf"))


tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg))

# plfl1 = drawPosesLandms(fg, spscale=1.0)
fg5a = deepcopy(fg)
tree5 = deepcopy(tree)


# drive second half

generateCanonicalFG_Circle(SIZE, fg=fg, biasTurn=-0.05, loopClosure=true, kappaOdo=3)

ensureAllInitialized!(fg)

# keep copy for later
fg_ = deepcopy(fg)



### SHOW Incremental Recycling--------------------------------------------------------


tree, smt, hist = solveTree!(fg, tree5, recordcliqs=ls(fg))



plotCirc10BA(fg_, fg, filepath=joinLogPath(fg, "circ$(SIZE)IncrBA.pdf"))


saveDFG(fg_, joinLogPath(fg,"fg_before_2nd"))
saveDFG(fg, joinLogPath(fg,"fg_after_2nd"))




## SHOW MARGINALIZATION---------------------------------------------------------


fg = deepcopy(fg_)

defaultFixedLagOnTree!(fg, round(Int, SIZE/2))

# use slightly special variable ordering
vo = getEliminationOrder(fg)
filter!(x->!(x in [:x5;:x4]), vo)
push!(vo, :x4)
push!(vo, :x5)

tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg), variableOrder=vo)

## Plot ellipses to illustrate covariance fit

plotCirc10BA(fg_, fg, filepath=joinLogPath(fg, "circ$(SIZE)MargBA.pdf"))


saveDFG(fg_, joinLogPath(fg,"fg_before_2nd"))
saveDFG(fg, joinLogPath(fg,"fg_after_2nd"))




## SHOW BOTH MARG AND INCR TOGETHER---------------------------------------------


fg = deepcopy(fg_)
defaultFixedLagOnTree!(fg, round(Int, SIZE/2))


tree, smt, hist = solveTree!(fg, tree5, recordcliqs=ls(fg))


plotCirc10BA(fg_, fg, filepath=joinLogPath(fg, "circ$(SIZE)ReclBA.pdf"))








## draw arcs

plotCirc10BA(fg5a, fg5a, filepath=joinLogPath(fg, "circ5A.pdf"))



tfg = initfg()

# add a starting point (skipping prior for brevity)
addVariable!(tfg, :a, Pose2)
manualinit!(tfg, :a, 0.01*randn(3,100))

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
