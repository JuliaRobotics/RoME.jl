
# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME

#  Do some plotting
using RoMEPlotting


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
fg.solverParams.isfixedlag = true
fg.solverParams.qfl = 15
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


getSolverParams(fg).drawtree = true

# debugging options
# getSolverParams(fg).multiproc = false
# getSolverParams(fg).async = true
# getSolverParams(fg).dbg = true
# getSolverParams(fg).limititers = 50
# getSolverParams(fg).downsolve = false


tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg))



pl = plotBeehive_6(fg, meanmax=:mean)
pl = plotBeehive_6(fg, meanmax=:max)

# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); @async run(`eog /tmp/test.svg`)

## fun additions
# fetchCliqTaskHistoryAll!(smt, hist)
# assignTreeHistory!(tree, hist)
# printCliqHistorySummary(hist, tree, :x1)
# makeCsmMovie(fg, tree)

# drawGraphCliq(hist, )

# plotTreeProductDown(fg, tree, :x4, :x4, levels=2)

# ## TEMP DEV CODE
#
# using Gadfly, Fontconfig, Cairo
# eo = getEliminationOrder(fg)
# tree = buildTreeFromOrdering!(fg, eo)
# drawTree(tree, show=true, imgs=true)
# writeGraphPdf(fg, show=true)
#
# dwinmsgs = prepCliqInitMsgsDown!(fg, tree, getCliq(tree, :x1), getCliq(tree,:x2), dbgnew=true)

# ## TEMP END


## hex 2

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br, autoinit=false )


posecount = driveHex(fg, posecount)


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )

# writeGraphPdf(fg,show=true)


# getSolverParams(fg).downsolve = false



tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)


# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)
pl = plotBeehive_6(fg)

fgc = deepcopy(fg)

deleteVariable!(fgc,:l2)

plotBeehive_6(fgc)

# new sighting on :l0

addVariable!(fg, :l0, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x5; :l0], p2br, autoinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x12; :l0], p2br2, autoinit=false )




# tree, smt, hist = solveTree!(fg, tree, recordcliqs=ls(fg))



## hex 3

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l3, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l3], p2br, autoinit=false )



posecount = driveHex(fg, posecount)



# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l3], p2br2, autoinit=false )


p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x19; :l0], p2br2, autoinit=false )

# writeGraphPdf(fg)


# getSolverParams(fg).dbg = true

tree, smt, hist = solveTree!(fg, tree, recordcliqs=ls(fg))




pl = plotBeehive_6(fg)
# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)

#
plotLocalProduct(fg, :l2)

plotTreeProductUp(fg,tree,:l2)

drawTree(tree, imgs=true)
drawTree(tree, imgs=false)

plotCliqUpMsgs(fg, tree,:l2)

ls(fg, :l2)
#
# pts = approxConv(fg, :x12l0f1, :l0)
#
# plotKDE(kde!(pts))
# plotPose(fg, :x12)
#
# plotKDE(fg, [:x12, :l0], dims=[1;2])


## hex 4

posecount = offsetHexLeg(fg, posecount, direction=:right)

posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)

tree, smt, hist = solveTree!(fg, tree)



drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)





## hex 5

posecount = offsetHexLeg(fg, posecount, direction=:right)

posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)


tree, smt, hist = solveTree!(fg, tree)






## hex 6

posecount = offsetHexLeg(fg, posecount, direction=:right)

posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)






## hex 7

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)

tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)



resetVariableAllInitializations!(fg)


unmarginalizeVariablesAll!(fg)


tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)





## hex 8

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)

tree, smt, hist = solveTree!(fg, tree)






## hex 9

posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


# writeGraphPdf(fg)

tree, smt, hist = solveTree!(fg, tree)







## hex 10

posecount = offsetHexLeg(fg, posecount, direction=:left)
posecount = offsetHexLeg(fg, posecount, direction=:right)


posecount = driveHex(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l1], p2br2, autoinit=false )


writeGraphPdf(fg)

tree, smt, hist = solveTree!(fg, tree)







0






using RoMEPlotting


drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)





## dev auto init after hex 3


printCliqHistorySummary(tree, :x16)

cliq = getCliq(tree, :x16)
cid = cliq.index
ste = 15
sfg = hist[cid][ste][4].cliqSubFg
cliqh = hist[cid][ste][4].cliq

printCliqSummary(sfg, cliq)

pl = plotBeehive_6(sfg)


stuff = sandboxCliqResolveStep(tree, :x16, 15)


cond = getSolveCondition(cliqh)
notify(cond)

cond.waitq = Any[]

isready(cond)

wait(getSolveCondition(cliqh))

notify(getSolveCondition(cliq))
IIF.notifyCliqUpInitStatus!(cliq, :downsolved)

getCliqStatus(cliq)

writeGraphPdf(sfg, show=true)


pl = plotBeehive_6(fg, to=15)




# leaf cliq

printCliqHistorySummary(tree, :l0)

cliq = getCliq(tree, :l0)
cid = cliq.index
ste = 9
sfg = hist[cid][ste][4].cliqSubFg
cliqh = hist[cid][ste][4].cliq


printCliqSummary(sfg, cliqh)

writeGraphPdf(sfg, show=true)

plotBeehive_6(sfg)


ls(sfg, :x19)

pts = approxConv(sfg, :x19l0f1, :x19)

X19 = manikde!(pts, Pose2().manifolds)

plotPose(Pose2(), X19)


pts[3,:] = 2*pi*rand(100).-pi
setValKDE!(getVariable(sfg, :x19), X19)


plotTreeProductDown(fg,tree,:l2, :l2)


## debug down on hex1, x4

drawTree(tree)



printCliqHistorySummary(tree,:x4)



csmcStep19 = getCliqSolveHistory(tree,:x4)[19][4]

drawGraph(csmcStep19.cliqSubFg, show=true)

getSolverParams(csmcStep19.dfg).dbg = true

stuff = sandboxCliqResolveStep(tree,:x4,19)


plotKDE(map(x->getKDE(x, :x4),[fg; csmcStep19.cliqSubFg; stuff[4].dfg]), levels=2 )

drawPosesLandms(fg,meanmax=:mean)



plotKDE([getKDE(csmcStep19.cliqSubFg, :x3);getKDE(csmcStep19.cliqSubFg, :x4); getKDE(csmcStep19.cliqSubFg, :x5)], levels=1)


plotLocalProduct(csmcStep19.cliqSubFg)


# what are down messages into (:x4|)

plotCliqDownMsgs(tree,:x3, existing=pl, dims=[1;2]);


plotTreeProductDown(csmcStep19.cliqSubFg, csmcStep19.tree,:x4,:x5)




## what is going on with x3 and x5 in root


printCliqHistorySummary(tree, :x3)
csmc = getCliqSolveHistory(tree, :x3)[9][4]


getSolverParams(csmc.dfg).dbg = true

stuff = sandboxCliqResolveStep(tree,:x3,9)


plotTreeUpMsgs(fg, tree, :x5) #plotCliqUpMsgs

plotCliqDownMsgs(tree,:x3, existing=pl, dims=[1;2]);


plotTreeProductUp(fg, tree, :x5)

#
