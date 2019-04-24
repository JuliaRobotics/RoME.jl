
# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME

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



## Manually organized a non-async initialization sequence

docliqs = [
5;
3;
4;
6;
2;
5;
1
]



# resetTreeCliquesForUpSolve!(tree)



# fg_bu = deepcopy(fg)
# tree_bu = deepcopy(tree)
# drawTree(tree_bu)


Niters = 100




ts5 = @async begin
cliq = tree.cliques[5]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=Niters )
end

ts3 = @async begin

cliq = tree.cliques[3]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=Niters )

end

ts4 = @async begin
cliq = tree.cliques[4]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=Niters )
end

ts6 = @async begin
cliq = tree.cliques[6]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=Niters )
end

ts2 = @async begin
cliq = tree.cliques[2]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=Niters )
end


ts1 = @async begin
cliq = tree.cliques[1]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=Niters )
end




ts1
ts2
ts3
ts4
ts5
ts6



drawTree(tree)


##








at = initInferTreeUp!(fg, tree, drawtree=true)




ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
# downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=true);
downMsgPassingRecursive(ett,N=100, dbg=false, drawpdf=true);







# Initialize :l1 numerical values but do not rerun solver
# ensureAllInitialized!(fg)
# pl = drawPosesLandms(fg, meanmax=:max)
# Gadfly.draw(Gadfly.PDF("/tmp/test2.pdf"),pl)  # or PNG(...)
# @async run(`evince /tmp/test2.pdf`)
