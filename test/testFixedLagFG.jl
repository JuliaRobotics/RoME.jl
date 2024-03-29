# example for fixed lag operation

using Test
using RoME
# using LinearAlgebra
# , Distributions


# TODO -- use this for cleanup
function driveSomeMore!(fg::AbstractDFG, idx)
    for i in idx:(idx+5)
      psym = Symbol("x$i")
      nsym = Symbol("x$(i+1)")
      addVariable!(fg, nsym, Pose2)
      pp = Pose2Pose2(MvNormal([10.0;0;pi/3], diagm([0.1;0.1;0.1].^2)))
      addFactor!(fg, [psym;nsym], pp )
    end
    return Symbol("x$(idx+5)")
end


##

@testset "test basic fixed lag operations..." begin
##

# start with an empty factor graph object
fg = initfg()
# IIF.getSolverParams(fg).drawtree = true
# IIF.getSolverParams(fg).showtree = true

# Set up a quasi fixed-lag horizon of 8 nodes and enable the fixed-lag solving.
# If the graph grows over 8 nodes, the older nodes will be frozen to limit the computational window.
IIF.getSolverParams(fg).qfl = 6
IIF.getSolverParams(fg).isfixedlag = true # limit up solve
IIF.getSolverParams(fg).limitfixeddown = true # also limit down solve


## 1. Drive around in a hexagon
# Add the first pose :x0
println("STEP 1: Driving around a bit")
addVariable!(fg, :x0, Pose2)
# Add at a fixed location PriorPose2 to pin :x0 to a starting location
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I,3,3))) )
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addVariable!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  addFactor!(fg, [psym;nsym], pp )
end

# Add node linking initial pose with a bearing range measurement landmark
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br)

## 2. Solve graph when shorter than fixed length - should solve full session.
println("STEP 2: Solve graph when shorter than fixed length")

# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true
tree = solveTree!(fg)

# Add another node when it comes around again, linking the node with the initial landmark
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br)


# 3. Drive a couple more, longer than fixed lag window
println("STEP 3: Drive a couple more, longer than fixed lag window")
for i in 6:11
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addVariable!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  addFactor!(fg, [psym;nsym], pp )
end

## At this point our window is 8 nodes, but our graph consists of 13 nodes.
## Next, freezing nodes beyond our fixed-lag horizon.

# Back up the graph
fgOriginal = deepcopy(fg)
getSolverParams(fgOriginal).isfixedlag = false

# Back up data from these two poses so we can compare them once we solve again.
X5 = deepcopy(getVal(fg, :x5))
X7 = deepcopy(getVal(fg, :x7))

@test isInitialized(fg, :x5)

fifoFreeze!(fg)
@test isMarginalized(fg, :x5)

# getSolverParams(fg).drawtree=true
# getSolverParams(fg).showtree=true
# getSolverParams(fg).dbg=true
# getSolverParams(fg).async=true

# Now solve again, which will freeze vertices < 5
println("STEP 4: Solve graph when shorter than fixed length, and show time to solve")
tree = solveTree!(fg, tree, verbose=true); # , recordcliqs=ls(fg));


# Confirm that the initial nodes (x0 - x5) are frozen.
@test isMarginalized(fg, :x5)
@test !isMarginalized(fg, :x7)

# X5 should be exactly same
# X6 should be different
X5cmp = deepcopy(getVal(fg, :x5))
X7cmp = deepcopy(getVal(fg, :x7))
for i in 1:length(X5)
  @test isapprox(submanifold_component(X5[i],1), submanifold_component(X5cmp[i],1))  #Frozen
  @test isapprox(submanifold_component(X5[i],2), submanifold_component(X5cmp[i],2))  #Frozen
end
for i in 1:length(X7)
  @test !isapprox(submanifold_component(X7[i],1), submanifold_component(X7cmp[i],1)) #Recalculated
  @test !isapprox(submanifold_component(X7[i],2), submanifold_component(X7cmp[i],2)) #Recalculated
end

##
end


#
