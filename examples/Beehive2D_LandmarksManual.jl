
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
fg.qfl = 25
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


# # writeGraphPdf(fg, show=true)
# tree = wipeBuildNewTree!(fg)
#
# cliq = tree.cliques[4]
# ids = IncrementalInference.getCliqFrontalVarIds(cliq)
# @show syms = map(d->getSym(fg, d), ids)
#


tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)


tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, downsolve=false , recordcliqs=[:x12;] )
tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true, upsolve=false ) #, recordcliqs=[:x3;]



# cliq = whichCliq(tree, :x3)
# hist = getData(cliq).statehistory

drawTree(tree)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); @async run(`eog /tmp/test.svg`)




## hex 2

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br, autoinit=false )


posecount = driveHex(fg, posecount, steps=5)

#
# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )

# writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg, drawpdf=true)



tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)  #, recordcliqs=[:x12;], limititers=30  , downsolve=false)


# prevx12hist = deepcopy(getCliqSolveHistory(tree, :x12))



# getData(whichCliq(tree, :x13)).statehistory = Vector{Tuple{Int, Function, CliqStateMachineContainer}}()





initInferTreeUp!(fg, tree, drawtree=true, limititers=20, recordcliqs=[:x12;])







drawTree(tree)

cliq = whichCliq(tree, :x12)
getCliqStatus(cliq)
hist = getCliqSolveHistory(tree, :x12)












# attempt redo of init up solve
iter = 20
csmc = deepcopy(hist[iter][3])


hist[iter][2](csmc)



prnt = getParent(csmc.tree, csmc.cliq)[1]
dwinmsgs = prepCliqInitMsgsDown!(csmc.cliqSubFg, csmc.tree, prnt)




cliq = whichCliq(tree, :x12)
prnt = getParent(tree, cliq)[1]
dwinmsgs = prepCliqInitMsgsDown!(fg, tree, prnt)











cliq = whichCliq(tree, :x12)
cliqInitSolveUpByStateMachine!(fg, tree, cliq, drawtree=true, limititers=20, recordhistory=true )



cliq = whichCliq(tree, :x13)
getCliqStatus(cliq)
cliqInitSolveUpByStateMachine!(fg, tree, cliq, drawtree=true, limititers=20, recordhistory=true )




cliq = whichCliq(tree,:x12)
prnt = getParent(tree, cliq)[1]
dwinmsgs = prepCliqInitMsgsDown!(fg, tree, prnt)






# import IncrementalInference: attemptCliqInitUp_StateMachine, doCliqAutoInitUp!, notifyCliqUpInitStatus!, initInferTreeUp!

function initInferTreeUp!(fgl::FactorGraph,
                          treel::BayesTree;
                          drawtree::Bool=false,
                          N::Int=100,
                          limititers::Int=-1,
                          recordcliqs::Vector{Symbol}=Symbol[] )
  #
  # revert :downsolved status to :initialized in preparation for new upsolve
  resetTreeCliquesForUpSolve!(treel)
  setTreeCliquesMarginalized!(fgl, treel)
  drawtree ? drawTree(treel, show=false) : nothing

  # queue all the tasks
  alltasks = Vector{Task}(undef, length(treel.cliques))
  cliqHistories = Dict{Int,Vector{Tuple{Int, Function, CliqStateMachineContainer}}}()
  @sync begin
    if !isTreeSolved(treel, skipinitialized=true)
      # duplicate int i into async (important for concurrency)
      for i in 1:length(treel.cliques)
        alltasks[i] = @async begin
          clst = :na
          cliq = treel.cliques[i]
          ids = getCliqFrontalVarIds(cliq)
          syms = map(d->getSym(fgl, d), ids)
          recordthiscliq = length(intersect(recordcliqs,syms)) > 0
          try
            history = cliqInitSolveUpByStateMachine!(fgl, treel, cliq, drawtree=drawtree, limititers=limititers, recordhistory=recordthiscliq )
            cliqHistories[i] = history
            clst = getCliqStatus(cliq)
            # clst = cliqInitSolveUp!(fgl, treel, cliq, drawtree=drawtree, limititers=limititers )
          catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
            error(err)
          end
          # if !(clst in [:upsolved; :downsolved; :marginalized])
          #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
          # end
        end # async
      end # for
    end # if
  end # sync

  # post-hoc store possible state machine history in clique (without recursively saving earlier history inside state history)
  for i in 1:length(treel.cliques)
    if haskey(cliqHistories, i)
      getData(treel.cliques[i]).statehistory=cliqHistories[i]
    end
  end

  return alltasks
end







drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)





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

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)


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





## hex 4

posecount = offsetHexLeg(fg, posecount, direction=:right)


# Add landmarks with Bearing range measurements
addVariable!(fg, :l4, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l4], p2br, autoinit=false )



posecount = driveHex(fg, posecount)


writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)




# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l4], p2br2, autoinit=false )




## Special sighting


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x19; :l0], p2br2, autoinit=false )



p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x26; :l0], p2br2, autoinit=false )






## hex 5

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l5, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l5], p2br, autoinit=false )


posecount = driveHex(fg, posecount)



writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") #|| @async run(`eog /tmp/test.svg`)




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




writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)

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



writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)

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
