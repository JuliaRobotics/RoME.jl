# example for dead-reckoning branch on factor graph

# using Revise

using DistributedFactorGraphs
using IncrementalInference
using RoME

using Test


@testset "testing against solvable..." begin

# start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) ) )

# Drive around in a hexagon
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp )
end

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, tags=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br )


## async solving with dead-reckon branch

addVariable!(fg, :deadreckon_x0, Pose2, solvable=0)

drec = MutablePose2Pose2Gaussian(MvNormal(zeros(3), Matrix{Float64}(LinearAlgebra.I, 3,3)))

addFactor!(fg, [:x0; :deadreckon_x0], drec, solvable=0)

#
@test length(map( x->x.label, getVariables(fg, solvable=1))) == 8
@test length(map( x->x.label, getVariables(fg, solvable=0))) == 9
#
# # make sure
@test length(getEliminationOrder(fg, solvable=1)) == 8
# check default
@test length(getEliminationOrder(fg)) == 8

# default check
vo = getEliminationOrder(fg)
@test length(vo) == 8


tree = buildTreeFromOrdering!(fg,vo)

# syms = map( x->x.label, getVariables(fg, solvable=1))

# union(map(x->getCliqAllVarIds(x), map(n->getClique(tree, n), syms))...)


# getSolverParams(fg).dbg = true
# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true


# breaks because internal saveDFG somewhere is failing on bad converter for MutablePose2Pose2Gaussian
tree2 = solveTree!(fg, recordcliqs=ls(fg));

@test !isInitialized(fg, :deadreckon_x0)


val = accumulateFactorMeans(fg, [:x0deadreckon_x0f1])

@test norm(val - calcVariablePPE(fg, :x0).suggested) < 1e-3


end









#
