using Revise
using Neo4j # So that DFG initializes the database driver.
using RoME
using DistributedFactorGraphs
using Test

# start with an empty factor graph object
# fg = initfg()
cloudFg = CloudGraphsDFG{SolverParams}("localhost", 7474, "neo4j", "test",
    "testUser", "testRobot", "testSession",
    nothing,
    nothing,
    IncrementalInference.decodePackedType,
    IncrementalInference.rebuildFactorMetadata!,
    solverParams=SolverParams())
# cloudFg = GraphsDFG{SolverParams}(params=SolverParams())
# cloudFg = GraphsDFG{SolverParams}(params=SolverParams())
clearSession!!(cloudFg)
# cloudFg = initfg()

# Add the first pose :x0
x0 = addVariable!(cloudFg, :x0, Pose2)
IncrementalInference.compareVariable(x0, getVariable(cloudFg, :x0))

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
prior = addFactor!(cloudFg, [:x0], PriorPose2( MvNormal([10; 10; 1.0/8.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) ) )
# retPrior = getFactor(cloudFg, :x0f1)
# Do the check
# IncrementalInference.compareFactor(prior, retPrior)
# Testing

# retPrior.data.fnc.cpt = prior.data.fnc.cpt
# # This one
# prior.data.fnc.cpt[1].factormetadata
# deserialized: Any[Pose2(3, String[], (:Euclid, :Euclid, :Circular))]
# vs.
# original: Pose2[Pose2(3, String[], (:Euclid, :Euclid, :Circular))]

# Drive around in a hexagon in the cloud
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(cloudFg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(cloudFg, [psym;nsym], pp )
end

# Right, let's copy it into local memory for solving...
localFg = GraphsDFG{SolverParams}(params=SolverParams())
DistributedFactorGraphs._copyIntoGraph!(cloudFg, localFg, union(getVariableIds(cloudFg), getFactorIds(cloudFg)), true)
# Some checks
@test symdiff(getVariableIds(localFg), getVariableIds(cloudFg)) == []
@test symdiff(getFactorIds(localFg), getFactorIds(cloudFg)) == []
@test isFullyConnected(localFg)
# Show it
toDotFile(localFg, "/tmp/localfg.dot")

# Alrighty! At this point, we should be able to solve locally...
# perform inference, and remember first runs are slower owing to Julia's just-in-time compiling
# Can do with graph too!
tree, smt, hist = solveTree!(localFg)

wipeBuildNewTree!(localFg)
tree, smt, hist = solveTree!(localFg, tree) # Recycle
# batchSolve!(localFg, drawpdf=true, show=true)
# Erm, whut? Error = mcmcIterationIDs -- unaccounted variables

# Trying new method.
tree, smtasks = batchSolve!(localFg, treeinit=true, drawpdf=true, show=true,
                            returntasks=true, limititers=50,
                            upsolve=true, downsolve=true  )

#### WIP and general debugging

# Testing with GenericMarginal
# This will not work because GenericMarginal *shouldn't* really be persisted.
# That would mean we're decomposing the cloud graph...
# genmarg = GenericMarginal()
# Xi = [getVariable(fg, :x0)]
# addFactor!(fg, Xi, genmarg, autoinit=false)

# For Juno/Jupyter style use
pl = drawPoses(localFg, meanmax=:mean)
plotPose(fg, :x6)
# For scripting use-cases you can export the image
Gadfly.draw(Gadfly.PDF("/tmp/test1.pdf", 20cm, 10cm),pl)  # or PNG(...)


# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br )


# Initialize :l1 numerical values but do not rerun solver
ensureAllInitialized!(fg)
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test2.pdf", 20cm, 10cm),pl)  # or PNG(...)


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br2 )


# solve
batchSolve!(fg, drawpdf=true)


# redraw
pl = drawPosesLandms(fg, meanmax=:mean)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)




#
