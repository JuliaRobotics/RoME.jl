
# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME
# import IncrementalInference: CliqStateMachineContainer

#  Do some plotting
# using RoMEPlotting


## start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([0.0; 0.0; 0.0],
                                           Matrix(Diagonal([0.1;0.1;0.05].^2))) ), autoinit=false )

# Drive around in a hexagon
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addVariable!(fg, nsym, Pose2)
  pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
  addFactor!(fg, [psym;nsym], pp, autoinit=false )
end

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false )


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x6; :l1], p2br2, autoinit=false )




tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=false)






cliq = whichCliq(tree, :x0)
history6 = cliqInitSolveUpByStateMachine!(fg, tree, cliq,
      drawtree=true, limititers=20, recordhistory=true)


drawTree(tree)




cliq = whichCliq(tree, :x2)
history5 = cliqInitSolveUpByStateMachine!(fg, tree, cliq,
      drawtree=true, limititers=30, recordhistory=true)

drawTree(tree)



cliq = whichCliq(tree, :x1)
history2 = cliqInitSolveUpByStateMachine!(fg, tree, cliq,
      drawtree=true, limititers=30, recordhistory=true)

drawTree(tree)


## init cliq 4

cliq = whichCliq(tree, :x4)

# csmc = CliqStateMachineContainer(fg, tree, cliq, initfg(), true, false, false, true, true)
csmc = CliqStateMachineContainer(fg, initfg(), tree, cliq, getParent(tree, cliq), getChildren(tree, cliq), true, false, false, true, true)

statemachine = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)






cliq3 = whichCliq(tree, :x6)
# csmc3 = CliqStateMachineContainer(fg, tree, cliq3, initfg(), true, false, false, true, true)
csmc3 = CliqStateMachineContainer(fg, initfg(), tree, cliq3, getParent(tree, cliq3), getChildren(tree, cliq3), true, false, false, true, true)


statemachine3 = StateMachine{CliqStateMachineContainer}(next=isCliqUpSolved_StateMachine)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)

# statemachine3.history







# statemachine.next
# blockCliqUntilChildrenHaveUpStatus

statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)
statemachine(csmc, verbose=true, recordhistory=true)

drawTree(tree)

# statemachine.history







statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)
statemachine3(csmc3, verbose=true, recordhistory=true)



drawTree(tree)


cliq = whichCliq(tree, :x3)
history1 = cliqInitSolveUpByStateMachine!(fg, tree, cliq,
      drawtree=true, limititers=30, recordhistory=true)

drawTree(tree)





using RoMEPlotting


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); @async run(`eog /tmp/test.svg`)



0
##



ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=true);


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg")






#
