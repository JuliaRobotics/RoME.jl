
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

# writeGraphPdf(fg, show=true)

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=false)


## Manually do tree based initialization


xx, ll = ls(fg)
for var in union(xx,ll)
  @show var, isInitialized(fg, var)
end

# fgc = deepcopy(fg)





## cliq 6

cliq = tree.cliques[6]
# cliq.attributes["label"]
# getCliqInitVarOrderUp(cliq)
syms = Symbol[getSym(fg, varid) for varid in getCliqAllVarIds(cliq)]

sfg = buildSubgraphFromLabels(fg, syms)
# writeGraphPdf(sfg, show=true)

doCliqAutoInitUp!(sfg, tree, cliq)
frsyms = Symbol[getSym(sfg, varid) for varid in getCliqFrontalVarIds(cliq)]
transferUpdateSubGraph!(sfg, fg, frsyms)



# isready(getData(cliq).initUpChannel)
# isCliqInitialized(cliq)



## cliq 5

# should be initializing with downward marginal message on x1 and x3

cliq = tree.cliques[5]
cliq.attributes["label"]
# getCliqInitVarOrderUp(cliq)
syms = Symbol[getSym(fg, varid) for varid in getCliqAllVarIds(cliq)]

sfg = buildSubgraphFromLabels(fg, syms)
# writeGraphPdf(sfg, show=true)


doCliqAutoInitUp!(sfg, tree, cliq)  # initstatus =




## first level up in tree to cliq2


pcliq = parentCliq(tree, :x0)[1]
# need some kind of blocking call till all siblings say the same
stdict = blockCliqUntilChildrenHaveUpStatus(tree, pcliq)


clid = 5
clst = stdict[clid]
for (clid, clst) in stdict

if clst == :needdownmsg

# could use parent sfg
prepCliqInitMsgsDown!(fg, tree, pcliq)

getCliqInitDownMsgs(pcliq)

doCliqAutoInitDown!(sfg, tree, cliq, dwinmsg)

end

end



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




##

ensureAllInitialized!(fg)


fg.isfixedlag = true
fg.qfl = 99


getData(fg, :x1).ismargin = true


using RoMEPlotting, Gadfly

XX1 = deepcopy(getKDE(fg, :x1))

plotPose(Pose2(), [XX1;]) |> SVG("/tmp/test.svg")
@async run(`inkscape /tmp/test.svg`)


plotPose(Pose2(), [getKDE(fg, :x1);]) |> SVG("/tmp/test2.svg")
@async run(`inkscape /tmp/test2.svg`)



# solve
batchSolve!(fg, drawpdf=true)

# redraw
pl = drawPosesLandms(fg, meanmax=:mean)
Gadfly.draw(Gadfly.PDF("/tmp/test3.pdf", 20cm, 10cm),pl)  # or PNG(...)




## testing subgraph

sfg = subgraphFromVerts(fg, [:x0;:x1;:l1], neighbors=2)

writeGraphPdf(sfg, show=true)





fgs = deepcopy(fg)




batchSolve!(fgs)



plotPose(fg, [:x1;:x2])


plotPose(fg, [:x1;:x2])



getKDE(sfg, :x1)



plotKDE([getKDE(fg, :x2); getKDE(fgs, :x2)], dims=[1;2], c=["red";"green"], levels=2)



getSym(fg, getCliqAllVarIds(whichCliq(tree, :x1))[2])


tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=true)

treeProductUp(fg, tree, :x0, :x0)






## Init Process summary`
#
# 1. trigger inits on all child cliques.
# 2. wait for (take!) response from all initUpChannel.
# 3. propagate initdownmsgs to any child cliq that needs down msgs.
# 4. initialize from down init msg and determine if further downward init msgs are required
## Note fully up-solvable only possible if all children completed up-solve


#
