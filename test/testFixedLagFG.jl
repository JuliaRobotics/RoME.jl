# example for fixed lag operation

using Test
using RoME
using LinearAlgebra
# , Distributions


# TODO -- use this for cleanup
function driveSomeMore!(fg::FactorGraph, idx)
    for i in idx:(idx+5)
      psym = Symbol("x$i")
      nsym = Symbol("x$(i+1)")
      addNode!(fg, nsym, Pose2)
      pp = Pose2Pose2(MvNormal([10.0;0;pi/3], diagm([0.1;0.1;0.1].^2)))
      addFactor!(fg, [psym;nsym], pp )
    end
    return Symbol("x$(idx+5)")
end


# start with an empty factor graph object
global fg = initfg() # FUTURE: (quasifixedwindow=50, autosolve=true)

# Set up a quasi fixed-lag horizon of 8 nodes and enable the fixed-lag solving.
# If the graph grows over 8 nodes, the older nodes will be frozen to limit the computational window.
fg.qfl = 8
fg.isfixedlag = true

@testset "test basic fixed lag operations..." begin

## 1. Drive around in a hexagon
# Add the first pose :x0
println("STEP 1: Driving around a bit")
addNode!(fg, :x0, Pose2)
# Add at a fixed location PriorPose2 to pin :x0 to a starting location
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I,3,3))) )
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addNode!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  addFactor!(fg, [psym;nsym], pp )
end

# Add node linking initial pose with a bearing range measurement landmark
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br)

## 2. Solve graph when shorter than fixed length - should solve full session.
println("STEP 2: Solve graph when shorter than fixed length")
IIF.batchSolve!(fg)

# 3. Drive a couple more, longer than fixed lag window
println("STEP 3: Drive a couple more, longer than fixed lag window")
for i in 6:11
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addNode!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  addFactor!(fg, [psym;nsym], pp )
end

# Add another node when it comes around again, linking the node with the initial landmark
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br)

## At this point our window is 8 nodes, but our graph consists of 13 nodes.
## Next, freezing nodes beyond our fixed-lag horizon.

# Back up the graph
fgOriginal = deepcopy(fg)
fgOriginal.isfixedlag = false

# Back up data from these two poses so we can compare them once we solve again.
X5 = deepcopy(getVal(fg, :x5))
X6 = deepcopy(getVal(fg, :x6))

# Now solve again, which will freeze vertices < 5
println("STEP 4: Solve graph when shorter than fixed length, and show time to solve")
@time IIF.batchSolve!(fg)
# fg.isfixedlag
# tofreeze = fg.fifo[1:(end-fg.qfl)]
# @test length(tofreeze) > 0
# IIF.setfreeze!.(fg, tofreeze)
#
# fifoFreeze!(fg)
# global tree = IIF.wipeBuildNewTree!(fg)
# inferOverTreeR!(fg, tree)

# Confirm that the initial nodes (x0 - x5) are frozen.
@test getData(fg, :x5).ismargin
@test getData(fg, :x6).ismargin == false

# X5 should be exactly same
# X6 should be different
X5cmp = deepcopy(getVal(fg, :x5))
X6cmp = deepcopy(getVal(fg, :x6))
@test X5 == X5cmp #Frozen
@test X6 != X6cmp #Recalculated

# Solve original graph with the to get time comparison
@time IIF.batchSolve!(fgOriginal)

end


#
