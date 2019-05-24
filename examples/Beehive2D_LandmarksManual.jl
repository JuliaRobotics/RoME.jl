
# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME

#  Do some plotting
using RoMEPlotting


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


function offsetHexLeg(fgl::FactorGraph, posecount::Int; direction=:right)
    psym = Symbol("x$(posecount-1)")
    nsym = Symbol("x$(posecount)")
    posecount += 1
    addVariable!(fgl, nsym, Pose2)
    pp = nothing
    if direction == :right
        pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    elseif direction == :left
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    end
    addFactor!(fgl, [psym; nsym], pp, autoinit=false )
    return posecount
end



## start with an empty factor graph object
fg = initfg()
fg.isfixedlag = true
fg.qfl = 20
posecount = 0

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)
posecount += 1


# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([0.0; 0.0; 0.0],
                                           Matrix(Diagonal([0.1;0.1;0.05].^2))) ), autoinit=false )

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false )



## hex 1

posecount = driveHex(fg, posecount)

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )




tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true,
                    returntasks=true, limititers=200, recordcliqs=[:x2;:x1;:x3],
                    downsolve=true, upsolve=true )


#


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); @async run(`eog /tmp/test.svg`)

0


## debug


# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)
# stuff = inferOverTree!(fg,tree,downsolve=false,drawpdf=true,treeinit=true,limititers=200,recordcliqs=Symbol[]  )



# tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true, downsolve=false ) #, recordcliqs=[:x3;] )
# tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true, upsolve=false ) #, recordcliqs=[:x3;]

printCliqHistorySummary(tree, :x2)



# hist = getCliqSolveHistory(tree, :x2)
# animateStateMachineHistoryByTime(hist, frames=100, folder="animatestate")



animateCliqStateMachines(tree, [:x2;:x1;:x3], frames=100)

# drawStateMachineHistory(hist, show=true)




# stuff = sandboxCliqResolveStep(tree, :x2, 6)


0
##

# cliq = whichCliq(tree, :x3)
# hist = getData(cliq).statehistory

# drawTree(tree)






## hex 2

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br, autoinit=false )


posecount = driveHex(fg, posecount, steps=5)


## Adding here fails partial tree init issue
# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )


# drawCliqSubgraphUp(fg, tree, :x7)


##


tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true,
                  returntasks=true, limititers=100, recordcliqs=[:x8;:x5;:x11;:x12;:x13],
                  downsolve=false, upsolve=true )
                  #, skipcliqids=[1;] )


0

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");

drawTree(tree)



printCliqHistorySummary(tree,:x12)
printCliqHistorySummary(tree,:x13)
printCliqHistorySummary(tree,:x11)
#


0


#
#
# tree = wipeBuildNewTree!(fg, drawpdf=true, imgs=false)
# # drawTree(tree, imgs=false)
#
#
#
#
# skipids = [1;2;3;4;5;6;7]
# stuff, cliqhist = initInferTreeUp!(fg,tree;drawtree=true,limititers=200,recordcliqs=Symbol[:x8;:x5;:x11;:x12;:x13], skipcliqids=skipids )
#
#
#
# cliq7 = whichCliq(tree, :x9)
# csmc7 = CliqStateMachineContainer(fg, initfg(), tree, cliq7, getParent(tree, cliq7), getChildren(tree, cliq7), false, true, true)
# statemachine7 = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
# for i in 1:11
#   statemachine7(csmc7, verbose=true, recordhistory=true)
#   drawTree(tree)
# end
#
# ts7 = @async statemachine7(csmc7, verbose=true, recordhistory=true)
#
# # printCliqHistorySummary(statemachine7.history)
# # drawTree(tree)
#
#
#
#
# cliq3 = whichCliq(tree, :x8)
# csmc3 = CliqStateMachineContainer(fg, initfg(), tree, cliq3, getParent(tree, cliq3), getChildren(tree, cliq3), false, true, true)
# statemachine3 = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
# for i in 1:4
#   statemachine3(csmc3, verbose=true, recordhistory=true)
#   drawTree(tree)
# end
#
# ts3 = @async statemachine3(csmc3, verbose=true, recordhistory=true)
# # drawTree(tree)
#
# printCliqHistorySummary(statemachine3.history)
#
#
#
#
# cliq4 = whichCliq(tree, :l1)
# csmc4 = CliqStateMachineContainer(fg, initfg(), tree, cliq4, getParent(tree, cliq4), getChildren(tree, cliq4), false, true, true)
# statemachine4 = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
# for i in 1:10
#   statemachine4(csmc4, verbose=true, recordhistory=true)
#   drawTree(tree)
# end
#
#
#
# printCliqHistorySummary(statemachine4.history)
# # drawTree(tree)
#
#
#
# ## definitely in a loop here, while waiting for :x9 down init to complete
# for i in 1:13
#   statemachine3(csmc3, verbose=true, recordhistory=true)
#   drawTree(tree)
# end
#
#
# printCliqHistorySummary(statemachine3.history)
#
#
#
#
# cliq7 = whichCliq(tree, :x9)
# csmc7 = CliqStateMachineContainer(fg, initfg(), tree, cliq7, getParent(tree, cliq7), getChildren(tree, cliq7), false, true, true)
# statemachine7 = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
# for i in 1:12
#   statemachine7(csmc7, verbose=true, recordhistory=true)
#   drawTree(tree)
# end
#
#
# printCliqHistorySummary(statemachine7.history)
#
#
#
#
#
#
# cliq2 = whichCliq(tree, :x5)
# csmc2 = CliqStateMachineContainer(fg, initfg(), tree, cliq2, getParent(tree, cliq2), getChildren(tree, cliq2), false, true, true)
# statemachine2 = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
# for i in 1:4
#   statemachine2(csmc2, verbose=true, recordhistory=true)
#   drawTree(tree)
# end
#
# t2 = @async statemachine2(csmc2, verbose=true, recordhistory=true)
#
#
# ts3 = @async statemachine3(csmc3, verbose=true, recordhistory=true)
# ts3 = @async statemachine3(csmc3, verbose=true, recordhistory=true)
# ts3 = @async statemachine3(csmc3, verbose=true, recordhistory=true)
# ts3 = @async statemachine3(csmc3, verbose=true, recordhistory=true)
# ts3 = @async statemachine3(csmc3, verbose=true, recordhistory=true)
# printCliqHistorySummary(statemachine3.history)
# drawTree(tree)
#
#
# statemachine2(csmc2, verbose=true, recordhistory=true)
# statemachine2(csmc2, verbose=true, recordhistory=true)
# statemachine2(csmc2, verbose=true, recordhistory=true)
# statemachine2(csmc2, verbose=true, recordhistory=true)
# statemachine2(csmc2, verbose=true, recordhistory=true)
#
#
#
# ## what is going on with the state machine here??
# printCliqHistorySummary(statemachine2.history)
#
#
#
#
#
#
# animateCliqStateMachines(tree, [:x8;:x5;:x11;:x12;], frames=500)
#
#
# getCliqStatus(whichCliq(tree, :x8))
#


# tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, downsolve=false , recordcliqs=[:x12;] )
# tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, upsolve=false ) #, recordcliqs=[:x3;]
#
# tree = wipeBuildNewTree!(fg, drawpdf=true)
# at = initInferTreeUp!(fg, tree, drawtree=true, limititers=65, recordcliqs=[:x12;])








# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )



## hex 3

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l3, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l3], p2br, autoinit=false )



posecount = driveHex(fg, posecount)


# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)








# Add landmarks with Bearing range measurements

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l3], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x18; :l1], p2br2, autoinit=false )




# new sighting

addVariable!(fg, :l0, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x5; :l0], p2br, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x12; :l0], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x19; :l0], p2br2, autoinit=false )



tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg")  #|| @async run(`eog /tmp/test.svg`)





## end debugging




## hex 4

posecount = offsetHexLeg(fg, posecount, direction=:right)


# Add landmarks with Bearing range measurements
addVariable!(fg, :l4, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l4], p2br, autoinit=false )



posecount = driveHex(fg, posecount)


# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)

0






# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l4], p2br2, autoinit=false )




## Special sighting


# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [:x19; :l0], p2br2, autoinit=false )



p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x26; :l0], p2br2, autoinit=false )






## hex 5

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l5, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l5], p2br, autoinit=false )


posecount = driveHex(fg, posecount)



# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)
# drawTree(tree, imgs=true)



# ls(fg)
# ls(fg, :l5)

# deleteFactor!(fg, :x28l5f1)
# deleteVariable!(fg, :l5)


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l5], p2br2, autoinit=false )




## Add more loop closures signthings


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x4; :l5], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x25; :l2], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x32; :l3], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x33; :l0], p2br2, autoinit=false )




# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)





## hex 6

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l6, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l6], p2br, autoinit=false )


posecount = driveHex(fg, posecount)

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l6], p2br2, autoinit=false )


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x11; :l6], p2br2, autoinit=false )


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x39; :l4], p2br2, autoinit=false )



# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)





drawTree(tree, imgs=true)




## hex 7

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)






## hex 8

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)






## hex 9

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)







## hex 10

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)







0






using RoMEPlotting


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)





#
