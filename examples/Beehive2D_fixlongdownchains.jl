
# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME

#  Do some plotting
# using RoMEPlotting


function driveHex(fgl, posecount::Int)
    # Drive around in a hexagon
    for i in (posecount-1):(posecount-1+5)
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
fg.qfl = 15
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


writeGraphPdf(fg, show=true)



tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)





## hex 2

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)







## hex 3

posecount = offsetHexLeg(fg, posecount, direction=:right)

posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)


tree = wipeBuildNewTree!(fg, drawpdf=true)




## Manually organized a non-async initialization sequence


## First half of tree (easy half)

doorder = [:x10; :x8; :x6; :x0; :x4; :x2; :x9; :l1; :x1; :x3]
docliqs = map(x->whichCliq(tree, x).index, doorder)


for i in docliqs
cliq = tree.cliques[i]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )
if !(clst in [:upsolved; :downsolved; :marginalized])
  error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
end
end







## Second half of tree

# easy branch
doorder = [:x14; :x12; :x13]
docliqs = map(x->whichCliq(tree, x).index, doorder)

for i in docliqs
cliq = tree.cliques[i]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )
if !(clst in [:upsolved; :downsolved; :marginalized])
  error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
end
end




## the long down chain problem


doorder = [:x18; :x16; :x19; :x17]
docliqs = map(x->whichCliq(tree, x).index, doorder)


cliq = tree.cliques[docliqs[1]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )


cliq = tree.cliques[docliqs[2]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )


cliq = tree.cliques[docliqs[3]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

drawTree(tree)

cliq = tree.cliques[docliqs[4]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=2 )

drawTree(tree)

# blocking call
# cliq = tree.cliques[docliqs[1]]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

cliq = tree.cliques[docliqs[2]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

cliq = tree.cliques[docliqs[3]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

cliq = tree.cliques[docliqs[1]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

cliq = tree.cliques[docliqs[3]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

cliq = tree.cliques[docliqs[4]]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=2 )


cliq = tree.cliques[1]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=2 )





ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=true);




0






using RoMEPlotting


pl = drawPosesLandms(fg, meanmax=:mean)
Gadfly.draw(Gadfly.SVG("/tmp/test2.svg"),pl)  # or PNG(...)

@async run(`eog /tmp/test2.svg`)




#
