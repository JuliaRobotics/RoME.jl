
# using Revise
using Distributed
addprocs(3)

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME
@everywhere using RoME

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
# fg.isfixedlag = true
# fg.qfl = 20
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





tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true,
                            returntasks=true, limititers=20,
                            upsolve=true, downsolve=false )
0

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true,
                            returntasks=true, limititers=20,
                            upsolve=false, downsolve=true )
0

# @info "Do recursive down inference over tree"
# downMsgPassingRecursive(ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[]),
#                           N=100, dbg=false, drawpdf=true);



#

# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");  @async run(`eog /tmp/test.svg`)





## debug



# writeGraphPdf(fg, show=true)
tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)
# drawTree(tree, show=true)
# allsyms = getTreeAllFrontalSyms(fg, tree);  #allsyms = Symbol[]
#                           recordcliqs=allsyms,



# plotPose(fg, :x0);

# tx6 = @async solveCliqWithStateMachine!(fg,tree,:x6, recordhistory=true, iters=50)
# tx4 = @async solveCliqWithStateMachine!(fg,tree,:x4, recordhistory=true, iters=50)

tx2 = @async solveCliqWithStateMachine!(fg,tree,:x2, recordhistory=true, iters=50)
sm0, csmc0 = solveCliqWithStateMachine!(fg,tree,:x0, recordhistory=true, iters=50)


# hist = getCliqSolveHistory(tree,:x2)
stuff = fetch(tx2)
hist = stuff[1].history
printCliqHistorySummary(hist)



## Debug dev

import IncrementalInference: attemptCliqInitDown_StateMachine, infocsm, attemptCliqInitUp_StateMachine, prepCliqInitMsgsDown!


@assert hist[14][3] == attemptCliqInitDown_StateMachine
pre14 = hist[14][4]
writeGraphPdf(pre14.cliqSubFg, engine="neato")

new15 = sandboxStateMachineStep(hist,14)
# Juno.@enter sandboxStateMachineStep(hist,14)


prnt = getParent(tree, pre14.cliq)
# getData(prnt[1])

## HERE IS THE PROBLEM
dwinmsgs = prepCliqInitMsgsDown!(pre14.cliqSubFg, pre14.tree, prnt[1])

writeGraphPdf(pre14.cliqSubFg, show=true)




hasVariable(pre14.cliqSubFg, :x1)






# # make sure csmc0 agrees
# ## why is this guy empty??
# dwinmsgs = prepCliqInitMsgsDown!(csmc0.cliqSubFg, csmc0.tree, csmc0.tree.cliques[prnt[1].index])
# dwinmsgs = prepCliqInitMsgsDown!(csmc0.cliqSubFg, csmc0.tree, csmc0.tree.cliques[prnt[1].index])



# how did :x0 solve go

printCliqHistorySummary(sm0.history)

csmc0new7 = sandboxStateMachineStep(sm0.history,7)

getData(csmc0new7[4].tree.cliques[2])








0



function prepCliqInitMsgsDown!(fgl::G,
                               tree::BayesTree,
                               cliq::Graphs.ExVertex ) where G <: AbstractDFG
  #
  @info "$(current_task()) Clique $(cliq.index), prepCliqInitMsgsDown!"
  # get the current messages stored in the parent
  currmsgs = getCliqInitUpMsgs(cliq)
  @info "$(current_task()) Clique $(cliq.index), msg keys=$(collect(keys(currmsgs)))"

  # check if any msgs should be multiplied together for the same variable
  msgspervar = Dict{Symbol, Vector{BallTreeDensity}}()
  for (cliqid, msgs) in currmsgs
    @show cliqid, length(msgs)
    for (msgsym, msg) in msgs
      if !haskey(msgspervar, msgsym)
        msgspervar[msgsym] = Vector{BallTreeDensity}()
      end
      push!(msgspervar[msgsym], msg)
    end
  end

  @info "$(current_task()) Clique $(cliq.index), keys with msgs=$(collect(keys(msgspervar)))"

  # reference to default allocated dict location
  products = getData(cliq).downInitMsg
  # multiply multiple messages together
  for (msgsym, msgs) in msgspervar
    # check if this particular down message requires msgsym
    @show msgsym, DFG.hasVariable(fgl, msgsym)
    if DFG.hasVariable(fgl, msgsym) #haskey(fgl.IDs, msgsym)
      if length(msgspervar[msgsym]) > 1
        products[msgsym] = manifoldProduct(msgs, getManifolds(fgl, msgsym))
      else
        products[msgsym] = msgs[1]
      end
    else
      # not required, therefore remove from message to avoid confusion
      if haskey(products, msgsym)
        delete!(products, msgsym)
      end
    end
  end

  @info "$(current_task()) Clique $(cliq.index), product keys=$(collect(keys(products)))"
  return products
end

0


function attemptCliqInitDown_StateMachine(csmc::CliqStateMachineContainer)
  #
  # should never happen to
  setCliqDrawColor(csmc.cliq, "green")
  csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

  # initialize clique in downward direction
  # not if parent also needs downward init message
  infocsm(csmc, "8a, needs down message -- attempt down init")
  prnt = getParent(csmc.tree, csmc.cliq)[1]
  @show dwinmsgs = prepCliqInitMsgsDown!(csmc.cliqSubFg, csmc.tree, prnt)

  @show cliqst = doCliqInitDown!(csmc.cliqSubFg, csmc.cliq, dwinmsgs)
  # TODO: transfer values changed in the cliques should be transfered to the tree in proc 1 here.

  # TODO: maybe this should be here?
  setCliqStatus!(csmc.cliq, cliqst)

  # TODO move out
  children = getChildren(csmc.tree, csmc.cliq)
  if areCliqChildrenNeedDownMsg(children)
    # set messages if children :needdownmsg
    infocsm(csmc, "8a, doCliqInitDown! -- must set messages for future down init")
    # construct init's up msg to place in parent from initialized separator variables
    msg = prepCliqInitMsgsUp(csmc.cliqSubFg, csmc.cliq) # , tree,

    infocsm(csmc, "8a, putting fake upinitmsg in this cliq, msgs labels $(collect(keys(msg)))")
    # set fake up and notify down status
    setCliqUpInitMsgs!(csmc.cliq, csmc.cliq.index, msg)
    # setCliqStatus!(csmc.cliq, cliqst)
    setCliqDrawColor(csmc.cliq, "sienna")
    csmc.drawtree ? drawTree(csmc.tree, show=false) : nothing

    notifyCliqDownInitStatus!(csmc.cliq, cliqst)

    infocsm(csmc, "8a, after down init attempt, $cliqst.")
  end

  return attemptCliqInitUp_StateMachine
end






















using Graphs

drawStateMachineHistory(hist)


drawTree(tree, imgs=true)

cliq = whichCliq(tree, :x1)

cliqd = getData(cliq)



# fihs = filterHistAllToArray(tree, [:x0;:x1;:x2;:x3;:x4;:x6], slowCliqIfChildrenNotUpsolved_StateMachine)
# printCliqHistorySummary(fihs)


##

# drawTree(tree)

# printCliqHistorySummary(tree, :x2)
# cliq = whichCliq(tree, :x3)
# hist = getData(cliq).statehistory
# animateCliqStateMachines(tree, [:x2;:x1;:x3], frames=100)



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
# tree = wipeBuildNewTree!(fg)
# allsyms = getTreeAllFrontalSyms(fg, tree);  #allsyms = Symbol[]
#                 recordcliqs=allsyms,


tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true,
                  returntasks=true, limititers=50, recordcliqs=[:x5;:x11;:x12;:x13],
                  upsolve=true, downsolve=true  )
0

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg");

# drawTree(tree)



# fihs = filterHistAllToArray(tree, allsyms, slowCliqIfChildrenNotUpsolved_StateMachine)
# printCliqHistorySummary(fihs)





printCliqHistorySummary(tree,:x5)
printCliqHistorySummary(tree,:x11)
printCliqHistorySummary(tree,:x12)
printCliqHistorySummary(tree,:x13)
#
# # hist = getCliqSolveHistory(tree, :x12)
# # sandboxCliqResolveStep(tree, :x12, 1)
# animateCliqStateMachines(tree,[:x12;:x11;:x8],frames=100)
#
#
# 0
##







#
#
# tree = wipeBuildNewTree!(fg, drawpdf=true, imgs=false)
# # drawTree(tree, imgs=false)
#







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



# solve

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true) #, recordcliqs=allsyms)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)







## closures

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


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #; @async run(`eog /tmp/test.svg`)
# drawTree(tree, show=true)




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



p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x25; :l2], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x26; :l0], p2br2, autoinit=false )


tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)






## hex 5

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l5, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l5], p2br, autoinit=false )


posecount = driveHex(fg, posecount)



# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)
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
addFactor!(fg, [:x32; :l3], p2br2, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x33; :l0], p2br2, autoinit=false )




# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)





## hex 6

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l6, Point2, labels=["LANDMARK"])
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



# writeGraphPdf(fg)

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)


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



tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, returntasks=true)

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)



drawTree(tree, imgs=true)




## hex 7

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)

recordcliqs = Symbol[] #[:x29;:x44;:x38;:x47;:x49]

tree, smtasks = batchSolve!(fg, treeinit=true, drawpdf=true, show=true,
                            returntasks=true, recordcliqs=recordcliqs )
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



# writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)

drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); # @async run(`eog /tmp/test.svg`)




## add sightings
addVariable!(fg, :l8, Point2, labels=["LANDMARK"])
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
