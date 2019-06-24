
# using Revise
using Distributed
addprocs(6)

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME

@everywhere using RoME

#  Do some plotting
using RoMEPlotting

0

function driveHex(fgl, posecount::Int; steps::Int=5)
    # Drive around in a hexagon
    for i in (posecount-1):(posecount-1+steps)
        psym = Symbol("x$i")
        posecount += 1
        nsym = Symbol("x$(i+1)")
        addVariable!(fgl, nsym, Pose2)
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
        addFactor!(fgl, [psym;nsym], pp, autoinit=false )
    end

    return posecount
end


function offsetHexLeg(dfg::G, posecount::Int; direction=:right) where G <: AbstractDFG
    psym = Symbol("x$(posecount-1)")
    nsym = Symbol("x$(posecount)")
    posecount += 1
    addVariable!(dfg, nsym, Pose2)
    pp = nothing
    if direction == :right
        pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    elseif direction == :left
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    end
    addFactor!(dfg, [psym; nsym], pp, autoinit=false )
    return posecount
end



## start with an empty factor graph object
fg = initfg()

# fg.solverParams
# fg.solverParams.isfixedlag = true
# fg.solverParams.qfl = 20
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

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true



tree, smt, chi = solveTree!(fg, recordcliqs=[:x3;:x1])


# hist = getCliqSolveHistory(tree, :x1)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");  @async run(`eog /tmp/test.svg`)




## hex 2

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br, autoinit=false )


posecount = driveHex(fg, posecount, steps=5)



# fg.solverParams.async = true

tree, smt, chi = solveTree!(fg, tree)

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");  @async run(`eog /tmp/test.svg`)




cliq = whichCliq(tree, :x10)
# smt[9]


prnt = getParent(tree, cliq)[1]
dwinmsgs = prepCliqInitMsgsDown!(fg, tree, prnt)




# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )


## hex 3

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l3, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l3], p2br, autoinit=false )


posecount = driveHex(fg, posecount)






tree, smt, hist = solveTree!(fg, tree)


# solve
#
# tree2, smtasks = batchSolve!(fg, tree, incremental=true, dbg=true, drawpdf=true, show=true)
# tree = tree2

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");  @async run(`eog /tmp/test.svg`)







## closures

# Add landmarks with Bearing range measurements

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l3], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x18; :l1], p2br2, autoinit=false )




# new sighting

addVariable!(fg, :l0, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x5; :l0], p2br, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x12; :l0], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x19; :l0], p2br2, autoinit=false )





tree, smt, hist = solveTree!(fg, tree)


# tree2, smtasks = batchSolve!(fg, tree, dbg=true, drawpdf=true, show=true, incremental=true)
# tree=tree2


drawPosesLandms(fg, meanmax=:max) |> PDF("/tmp/test.pdf");  @async run(`evince /tmp/test.pdf`)
# drawTree(tree, show=true)





## hex 4

posecount = offsetHexLeg(fg, posecount, direction=:right)


# Add landmarks with Bearing range measurements
addVariable!(fg, :l4, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l4], p2br, autoinit=false )



posecount = driveHex(fg, posecount)


# writeGraphPdf(fg)


tree, smt, hist = solveTree!(fg, tree)

# tree2, smtasks = batchSolve!(fg, tree, dbg=true, drawpdf=true, show=true, incremental=true)
# tree=tree2




drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); @async run(`eog /tmp/test.svg`)


0



# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x27; :l4], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x25; :l2], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x26; :l0], p2br2, autoinit=false )




tree, smt, hist = solveTree!(fg, tree)



# tree2, smtasks = batchSolve!(fg, tree, incremental=true, dbg=true, drawpdf=true, show=true)
# tree = tree2


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");  @async run(`eog /tmp/test.svg`)



ls(fg, :l4)



## hex 5

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l5, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l5], p2br, autoinit=false )


posecount = driveHex(fg, posecount)





tree, smt, hist = solveTree!(fg, tree)



# tree2, smtasks = batchSolve!(fg, tree, incremental=true, dbg=true, drawpdf=true, show=true)
# tree = tree2


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)
# drawTree(tree, imgs=true)



# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l5], p2br2, autoinit=false )



## Add more loop closures signthings

# THIS IS WRONG
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [:x4; :l5], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x32; :l3], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x33; :l0], p2br2, autoinit=false )




tree, smt, hist = solveTree!(fg, tree)



# tree2, smtasks = batchSolve!(fg, tree, incremental=true, dbg=true, drawpdf=true, show=true)
# tree = tree2


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");  @async run(`eog /tmp/test.svg`)





## hex 6

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l6, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l6], p2br, autoinit=false )


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l6], p2br2, autoinit=false )
#
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [:x11; :l6], p2br2, autoinit=false )
#
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [:x39; :l4], p2br2, autoinit=false )




tree, smt, hist = solveTree!( fg, tree )


# tree2, smtasks = batchSolve!(fg, tree, incremental=true, dbg=true, drawpdf=true, show=true)
# tree = tree2


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)







# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x40; :l0], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l6], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x11; :l6], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x39; :l4], p2br2, autoinit=false )



tree2, smtasks = batchSolve!(fg, tree, incremental=true, dbg=true, drawpdf=true, show=true)
tree = tree2

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)|



drawTree(tree, imgs=true)




## hex 7

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)


# Add landmarks with Bearing range measurements
addVariable!(fg, :l7, Point2, labels=[:LANDMARK])
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x49; :l7], p2br2, autoinit=false )
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x2; :l7], p2br2, autoinit=false )



# writeGraphPdf(fg)

recordcliqs = [:x29;:x44;:x38;:x47;:x49]


tree2, smtasks = batchSolve!(fg, tree, incremental=false, dbg=true, drawpdf=true, show=true, recordcliqs=recordcliqs)
tree = tree2
0


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)




# hist = getCliqSolveHistory(tree, :x38)
# printCliqHistorySummary(tree, :x38)
# printCliqHistorySummary(tree, :x47)
#
# animateStateMachineHistoryByTime(hist)
#
#
#
# sandboxCliqResolveStep(tree,:x38,37)
# sandboxCliqResolveStep(tree,:x38,41)
#



## hex 8

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)



tree, smt, chi = solveTree!(fg, tree)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)




## add sightings
addVariable!(fg, :l8, Point2, labels=[:LANDMARK])
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x2; :l8], p2br2, autoinit=false )


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x49; :l8], p2br2, autoinit=false )



## more resightings

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x47; :l5], p2br2, autoinit=false )


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x48; :l6], p2br2, autoinit=false )


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x55; :l6], p2br2, autoinit=false )




tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)







# unmarginalizeVariablesAll!(fg)

plotPose(fg, :x22);

## hex 9

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )



tree, smt, chi = solveTree!(fg, tree)






## hex 10

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )








0






using RoMEPlotting


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)





#
