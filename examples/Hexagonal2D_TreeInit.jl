
using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
using Cairo, Fontconfig
using Gadfly

using RoME


## start with an empty factor graph object
fg = initfg()

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) ), autoinit=true )

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
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false )


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br2, autoinit=false )

##

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=true)


## Manually do tree based initialization


xx, ll = ls(fg)
for var in union(xx,ll)
  @show var, isInitialized(fg, var)
end

fgc = deepcopy(fg)

## cliq 6

cliq = tree.cliques[6]
cliq.attributes["label"]
getCliqInitVarOrderUp(cliq)


doCliqAutoInitUp!(fg, tree, cliq)
areCliqVariablesInitialized(fg, cliq)
getData(cliq).initialized
isCliqReadyInferenceUp(fg, tree, cliq)




## cliq 5

# should be initializing with downward marginal message on x1 and x3

cliq = tree.cliques[5]
cliq.attributes["label"]
getCliqInitVarOrderUp(cliq)


doCliqAutoInitUp!(fg, tree, cliq)



## cliq 2

cliq = tree.cliques[2]
cliq.attributes["label"]
getCliqInitVarOrderUp(cliq)


doCliqAutoInitUp!(fg, tree, cliq)


approxCliqMarginalUp!(fg,tree, x1, true)


# cliq 3

cliq = tree.cliques[3]
cliq.attributes["label"]
getCliqInitVarOrderUp(cliq)


doCliqAutoInitUp!(fg, tree, cliq)



# cliq 4

cliq = tree.cliques[4]
cliq.attributes["label"]
getCliqInitVarOrderUp(cliq)


doCliqAutoInitUp!(fg, tree, cliq)





##  Do some plotting
using RoMEPlotting


# Initialize :l1 numerical values but do not rerun solver
# ensureAllInitialized!(fg)
pl = drawPosesLandms(fg)
Gadfly.draw(Gadfly.PDF("/tmp/test2.pdf", 20cm, 10cm),pl)  # or PNG(...)




# solve
batchSolve!(fg, drawpdf=true)

# redraw
pl = drawPosesLandms(fg, meanmax=:mean)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)




## testing subgraph

sfg = subgraphFromVerts(fg, [:x0;:x1;:l1], neighbors=2)

writeGraphPdf(sfg, show=true)


#
