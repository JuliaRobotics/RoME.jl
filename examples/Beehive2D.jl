
using Distributed
addprocs(4)

# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using Dates
using RoME

import RoME: _driveHex!, _offsetHexLeg

@everywhere using RoME

#  Do some plotting
using RoMEPlotting
Gadfly.set_default_plot_size(35cm,25cm)

##

# function _driveHex!(fgl::AbstractDFG,
#                   posecount::Int;
#                   graphinit::Bool=false )
#     #

#     # Drive around in a hexagon
#     for i in (posecount):(posecount+5)
#         psym = Symbol("x$i")
#         posecount += 1
#         nsym = Symbol("x$(i+1)")
#         addVariable!(fgl, nsym, Pose2)
#         pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
#         addFactor!(fgl, [psym;nsym], pp, graphinit=graphinit )
#     end

#     return posecount+6
# end


# function _offsetHexLeg(fgl::G,
#                       posecount::Int; direction=:right,
#                       graphinit::Bool=false ) where G <: AbstractDFG
#     #
#     psym = Symbol("x$(posecount-1)")
#     nsym = Symbol("x$(posecount)")
#     posecount += 1
#     addVariable!(fgl, nsym, Pose2)
#     pp = nothing
#     if direction == :right
#         pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
#     elseif direction == :left
#         pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
#     end
#     addFactor!(fgl, [psym; nsym], pp, graphinit=graphinit )
#     return posecount
# end



## start with an empty factor graph object

fg = initfg()
getSolverParams(fg).isfixedlag = true
getSolverParams(fg).qfl = 20
posecount = 0

# Add the first pose :x0
addVariable!(fg, :x0, Pose2)
# posecount += 1


# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([0.0; 0.0; 0.0],
                                           Matrix(Diagonal([0.1;0.1;0.05].^2))) ), graphinit=false )

# Add landmarks with Bearing range measurements
addVariable!(fg, :l0, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x0; :l0], p2br, graphinit=false )




## hex 1

posecount = _driveHex!(fg, posecount)

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount)"); :l0], p2br2, graphinit=false )


# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x1; :l1], p2br, graphinit=false )



# draw figures for debugging
getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = false

# burn out USB drive instead of solid storage
# getSolverParams(fg).logpath="/media/dehann/temp/caesar/$(now())"

# debugging options
# getSolverParams(fg).multiproc = true
# getSolverParams(fg).async = true
# getSolverParams(fg).dbg = true
# getSolverParams(fg).limititers = 50
# getSolverParams(fg).downsolve = false



tree,_ , = solveTree!(fg)


pl = plotBeehive_6(fg, meanmax=:mean)
# pl = plotBeehive_6(fg, meanmax=:max)
# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg"); @async run(`eog /tmp/test.svg`)


## fun additions
# fetchCliqTaskHistoryAll!(smt, hist)
# printCliqHistorySummary(hist, tree[:l0])
# makeCsmMovie(fg, tree)



## hex 2

posecount = _offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount)"); :l2], p2br, graphinit=false )


posecount = _driveHex!(fg, posecount)


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount)"); :l2], p2br2, graphinit=false )

# drawGraph(fg,show=true)


# getSolverParams(fg).downsolve = false


tree = solveTree!(fg, tree)



# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)
pl = plotBeehive_6(fg, meanmax=:mean)



## new sighting on :l0

addVariable!(fg, :l0, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x5; :l0], p2br, graphinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x12; :l0], p2br2, graphinit=false )



# getSolverParams(fg).dbg = true

tree = solveTree!(fg, tree) #, recordcliqs=ls(fg))


pl = plotBeehive_6(fg, meanmax=:mean)

# getSolverParams(fg).dbg = false



## hex 3

posecount = _offsetHexLeg(fg, posecount, direction=:right, graphinit=true)


posecount = _driveHex!(fg, posecount, graphinit=true)



tree = solveTree!(fg, tree)


pl = plotBeehive_6(fg, meanmax=:mean)


##

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x19; :l0], p2br2, graphinit=false )

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x18; :l0], p2br2, graphinit=false )


# Add landmarks with Bearing range measurements
addVariable!(fg, :l3, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x14; :l3], p2br, graphinit=false )


p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x20; :l3], p2br, graphinit=false )



# getSolverParams(fg).dbg = true

tree = solveTree!(fg, tree)


pl = plotBeehive_6(fg, meanmax=:mean)
# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)





## hex 4

posecount = _offsetHexLeg(fg, posecount, direction=:right)

posecount = _driveHex!(fg, posecount)



# Add landmarks with Bearing range measurements
addVariable!(fg, :l4, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x21; :l4], p2br, graphinit=false )

p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x27; :l4], p2br, graphinit=false )



# drawGraph(fg)

tree = solveTree!(fg, tree)


pl = plotBeehive_6(fg, meanmax=:mean)


##

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x26; :l0], p2br2, graphinit=false )

# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x25; :l2], p2br2, graphinit=false )


tree = solveTree!(fg, tree, recordcliqs=ls(fg))


pl = plotBeehive_6(fg, meanmax=:mean)
# drawPosesLandms(fg, meanmax=:max) |> SVG("/tmp/test.svg") # || @async run(`eog /tmp/test.svg`)




## hex 5

posecount = _offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l5, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x28; :l5], p2br, graphinit=false )


p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x4; :l5], p2br, graphinit=false )


tree = solveTree!(fg, tree)

pl = plotBeehive_6(fg, meanmax=:mean)



## Weird initialization issues happen here


posecount = _driveHex!(fg, posecount, graphinit=true)

# # Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x34; :l5], p2br2, graphinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x32; :l3], p2br2, graphinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x33; :l0], p2br2, graphinit=false )



# drawGraph(fg)


tree = solveTree!(fg, tree)

pl = plotBeehive_6(fg, meanmax=:mean)




## hex 6

posecount = _offsetHexLeg(fg, posecount, direction=:right)

posecount = _driveHex!(fg, posecount, graphinit=true)

# l6 from x35, x41, x11
# Add landmarks with Bearing range measurements
addVariable!(fg, :l6, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x35; :l6], p2br, graphinit=false )

p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x41; :l6], p2br, graphinit=false )

p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x11; :l6], p2br, graphinit=false )



# # Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x40; :l0], p2br2, graphinit=false )

p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x39; :l4], p2br2, graphinit=false )


# drawGraph(fg)

tree = solveTree!(fg, tree)


pl = plotBeehive_6(fg, meanmax=:mean)




## complete solve

# dontMarginalizeVariablesAll!(fg)
#
# # getSolverParams(fg).dbg = false
# # getSolverParams(fg).multiproc = false
#
# tree = solveTree!(fg) # , recordcliqs=ls(fg)
#
# pl = plotBeehive_6(fg, meanmax=:mean)
#
#
# # plotLocalProduct(fg, :l6)
# getSolverParams(fg).isfixedlag = true
# getSolverParams(fg).qfl = 20



## hex 7

posecount = _offsetHexLeg(fg, posecount, direction=:left)
posecount = _offsetHexLeg(fg, posecount, direction=:right)


posecount = _driveHex!(fg, posecount, graphinit=true)


#
# # # Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x47; :l5], p2br2, graphinit=false )
#
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x48; :l6], p2br2, graphinit=false )



tree = solveTree!(fg, tree)


pl = plotBeehive_6(fg, meanmax=:mean)
# resetVariableAllInitializations!(fg)



## hex 8

posecount = _offsetHexLeg(fg, posecount, direction=:right)


posecount = _driveHex!(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l0], p2br2, graphinit=false )


# drawGraph(fg)

tree = solveTree!(fg, tree)


pl = plotBeehive_6(fg, meanmax=:mean)



## hex 9

posecount = _offsetHexLeg(fg, posecount, direction=:right)


posecount = _driveHex!(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l0], p2br2, graphinit=false )


# drawGraph(fg)

tree = solveTree!(fg, tree)



pl = plotBeehive_6(fg, meanmax=:mean)



## hex 10

posecount = _offsetHexLeg(fg, posecount, direction=:left)
posecount = _offsetHexLeg(fg, posecount, direction=:right)


posecount = _driveHex!(fg, posecount)

# # Add landmarks with Bearing range measurements
# p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
# addFactor!(fg, [Symbol("x$(posecount-1)"); :l0], p2br2, graphinit=false )


drawGraph(fg)

tree = solveTree!(fg, tree)




plotSLAM2D(fg, meanmax=:max) |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)







## additional dev code

# from hex1

getSolverParams(fg).upsolve = false
getSolverParams(fg).downsolve = true

stuff = solveCliq!( fg, tree, :x3 )
stuff = solveCliq!( fg, tree, :x1 )
stuff = solveCliq!( fg, tree, :x6 )
stuff = solveCliq!( fg, tree, :x4 )
stuff = solveCliq!( fg, tree, :x2 )
stuff = solveCliq!( fg, tree, :x0 )

fg2bd = loadDFG("/tmp/caesar/cliqSubFgs/cliq2/fg_beforedownsolve", Main)

drawGraph(fg2bd, show=true)





## dev auto init after hex 3


printCliqHistorySummary(tree, :x16)

cliq = getClique(tree, :x16)
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

drawGraph(sfg, show=true)


pl = plotBeehive_6(fg, to=15)




# leaf cliq

printCliqHistorySummary(tree, :l0)

cliq = getClique(tree, :l0)
cid = cliq.index
ste = 9
sfg = hist[cid][ste][4].cliqSubFg
cliqh = hist[cid][ste][4].cliq


printCliqSummary(sfg, cliqh)

drawGraph(sfg, show=true)

plotBeehive_6(sfg)


ls(sfg, :x19)

pts = approxConv(sfg, :x19l0f1, :x19)

X19 = manikde!(Pose2, pts)

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


plotKDE(map(x->getBelief(x, :x4),[fg; csmcStep19.cliqSubFg; stuff[4].dfg]), levels=2 )

drawPosesLandms(fg,meanmax=:mean)



plotKDE([getBelief(csmcStep19.cliqSubFg, :x3);getBelief(csmcStep19.cliqSubFg, :x4); getBelief(csmcStep19.cliqSubFg, :x5)], levels=1)


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



hist = getCliqSolveHistory(tree, :x6)

dfg = hist[11][4].dfg
subfg = hist[11][4].cliqSubFg
cliq = hist[11][4].cliq

printCliqHistorySummary(hist)




## DEBUG x19

plotKDE(fg, :x19, levels=3)

getClique(tree, :x19)

printCliqHistorySummary(tree, :x19)

sfg_ad = hist[4][12][4].cliqSubFg

drawGraph(sfg_ad, show=true)

plotLocalProduct(sfg_ad, :x19)

getVariableInferredDim(fg, :x19)
getVariableInferredDim(fg, :l0)


## END DEBUG




"""
    $SIGNATURES



Related

plotLocalProduct
"""
function plotLocalInverseProduct()

end


plotKDE(fg, :x41, levels=2)

plotLocalProduct(fg, :x41)






## DEBUG weird hex5

dontMarginalizeVariablesAll!(fg)
# getSolverParams(fg).isfixedlag = false
# getSolverParams(fg).qfl = 99999999

getSolverParams(fg).dbg = true
# getSolverParams(fg).multiproc = false

tree = solveTree!(fg, recordcliqs=ls(fg))

pl = plotBeehive_6(fg, meanmax=:mean)


plotLocalProduct(fg, :l5)

stuff = localProduct(fg,:l5)

setValKDE!(fg, :l5, stuff[1], false, stuff[5])

ls(fg, :x33)


saveDFG(fg, "/home/dehann/Documents/beehive_hex5_fail_9_13")


printCliqHistorySummary(tree,:l0)

l0bdfg = getCliqSolveHistory(tree,:l0)[11][4].dfg
l0adfg = getCliqSolveHistory(tree,:l0)[12][4].dfg

plotKDE(l0bdfg, :l0, levels=3)
plotKDE(l0adfg, :l0, levels=3)

getPoints(getBelief(l0bdfg, :l0))


getPoints(getBelief(l0adfg, :l0))

## DEBUG END






## DEBUG post 6

pl = plotKDE(fg, ls(fg, r"l"), levels=1)
pl = plotKDE(fg, ls(fg, r"x"), levels=1)

dontMarginalizeVariablesAll!(fg)
getSolverParams(fg).dbg=true
tree = solveTree!(fg, recordcliqs=ls(fg)) #, tree)

pl = plotBeehive_6(fg, meanmax=:mean)

cliq = getClique(tree, :l0)

printCliqHistorySummary(tree, :l0)

sfg31_ad = getCliqSolveHistory(tree, :l0)[12][4].cliqSubFg
fg31_11 = getCliqSolveHistory(tree, :l0)[12][4].dfg
fg31_12 = getCliqSolveHistory(tree, :l0)[12][4].dfg

plotKDE(sfg31_ad, :l0, levels=1)

plotLocalProduct(sfg31_ad, :l0)
plotLocalProduct(fg, :l0)

drawGraph(sfg31_ad)

plotKDE(sfg31_ad, sortDFG(ls(sfg31_ad)), levels=1)

plotKDE(fg31_11, sortDFG(ls(sfg31_ad)), levels=1)
plotKDE(fg31_12, sortDFG(ls(sfg31_ad)), levels=1)

plotCliqUpMsgs(fg, tree, :x5)

plotCliqDownMsgs(tree, :x5)


getDwnMsgs(tree, :x9)



ls(fg, :l6)

plotKDE(fg, :l6, levels=2)

plotLocalProduct(fg, :l6)

getSolverParams(fg).isfixedlag = true
getSolverParams(fg).qfl = 15

## DEBUG END



#
