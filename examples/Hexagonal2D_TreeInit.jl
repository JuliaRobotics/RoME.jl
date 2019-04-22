
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

tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=true)


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
writeGraphPdf(sfg, show=true)

doCliqAutoInitUp!(sfg, tree, cliq)
frsyms = Symbol[getSym(sfg, varid) for varid in getCliqFrontalVarIds(cliq)]
transferUpdateSubGraph!(fg, sfg, frsyms)



# isready(getData(cliq).initUpChannel)
# isCliqInitialized(cliq)



## cliq 5

# should be initializing with downward marginal message on x1 and x3


cliq = tree.cliques[5]
cliq.attributes["label"]
# getCliqInitVarOrderUp(cliq)
syms = getCliqAllVarSyms(fg, cliq)

sfg = buildSubgraphFromLabels(fg, syms)
writeGraphPdf(sfg, show=true)


doCliqAutoInitUp!(sfg, tree, cliq)  # initstatus =


# level up in tree to cliq2
pcliq = parentCliq(tree, :x0)[1]
# blocking call till all siblings have a usable status
stdict = blockCliqUntilChildrenHaveUpStatus(tree, pcliq)


clid = 5
clst = stdict[clid]


if clst == :needdownmsg

# initialize clique in downward direction
clst = doCliqInitDown!(sfg, tree, cliq)

end # :needdownmsg


if clst == :initialized

clst = doCliqAutoInitUp!(sfg, tree, cliq)

end

if clst == :upsolved

frsyms = Symbol[getSym(sfg, varid) for varid in getCliqFrontalVarIds(cliq)]
transferUpdateSubGraph!(fg, sfg, frsyms)

end # :initialized




IncrementalInference.drawTree(tree, show=true)







## cliq 2




clid = 2
cliq = tree.cliques[clid]
cliq.attributes["label"]
syms = getCliqAllVarSyms(fg, cliq)

# build a local subgraph for inference operations
sfg = buildSubgraphFromLabels(fg, syms)
writeGraphPdf(sfg, show=true)


tryonce = true

# upsolve delay loop
# while !areCliqChildrenAllUpSolved(tree, cliq) || tryonce

tryonce = false

# wait here until all children have a valid status
stdict = blockCliqUntilChildrenHaveUpStatus(tree, cliq)

proceed = true
# possible send down msg for pending child status
for (clid, clst) in stdict
  # :needdownmsg # 'send' downward init msg direction
  # :initialized # @warn "something might not be right with init of clid=$clid"
  clst != :upsolved ? (proceed = false) : nothing
end

# if all children are ready, proceed with this cliq initialization
if proceed
proceed = false

cliqst = getCliqStatus(cliq)

if cliqst == :needdownmsg
  # initialize clique in downward direction
  cliqst = doCliqInitDown!(sfg, tree, cliq)
end
if cliqst in [:initialized; :null]
  cliqst = doCliqAutoInitUp!(sfg, tree, cliq)
end
if cliqst == :upsolved
  frsyms = Symbol[getSym(sfg, varid) for varid in getCliqFrontalVarIds(cliq)]
  transferUpdateSubGraph!(fg, sfg, frsyms)
else
  @info "clique $(cliq.index) init waiting since it cannot fully up solve yet."
end

end


# end # while



drawTree(tree, show=true)




## Pause here to figure out what is happening with up messages





## continue with clique inference



# cliq 3

cliq = tree.cliques[3]
cliq.attributes["label"]


cliqInitSolveUp!(fg, tree, cliq )





# cliq 4

cliq = tree.cliques[4]
cliq.attributes["label"]


cliqInitSolveUp!(fg, tree, cliq )




drawTree(tree, show=true)






# cliq 4

cliq = tree.cliques[1]
cliq.attributes["label"]


cliqInitSolveUp!(fg, tree, cliq )



drawTree(tree, show=true)




## follow with downward inference

ett = ExploreTreeType(fg, tree, cliq, nothing, NBPMessage[])

IncrementalInference.downMsgPassingRecursive(ett, dbg=false, drawpdf=true);






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


using RoMEPlotting


plotPose(fg,:x2)
