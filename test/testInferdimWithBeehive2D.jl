
# using Revise

##

using RoME

#  Do some plotting
# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)


function driveHex(fgl::G, posecount::Int) where G <: AbstractDFG
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


function offsetHexLeg(fgl::G, posecount::Int; direction=:right) where G <: AbstractDFG
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
getSolverParams(fg).isfixedlag = true
getSolverParams(fg).qfl = 15
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

# writeGraphPdf(fg, show=true)


# debugging options
getSolverParams(fg).drawtree = true
getSolverParams(fg).multiproc = false
getSolverParams(fg).dbg = true
getSolverParams(fg).limititers = 50
# getSolverParams(fg).downsolve = false
getSolverParams(fg).async = true


# do inference over factor graph
tree, smt, hist = solveTree!(fg, recordcliqs=ls(fg))


fetchCliqTaskHistoryAll!(smt, hist)
assignTreeHistory!(tree, hist)






## test first cliq for inferdim

printCliqHistorySummary(tree, :x0)
# cliq6 = getCliq(tree, :x0)

sfg = hist[6][9][4].cliqSubFg

isInitialized(sfg,:x0)
getVariableInferredDim(sfg,:x0)
isInitialized(fg,:x0)
getVariableInferredDim(fg,:x0)

isInitialized(sfg,:x1)
getVariableInferredDim(sfg,:x1)
isInitialized(fg,:x1)
getVariableInferredDim(fg,:x1)

isInitialized(sfg,:l1)
getVariableInferredDim(sfg,:l1)
isInitialized(fg,:l1)
getVariableInferredDim(fg,:l1)


## check the up message sent to

sfg = hist[6][9][4].cliqSubFg
upmsg = prepCliqInitMsgsUp(sfg, getCliq(tree, :x0))
upMsg(tree, :x0)


dwnmsg = getDwnMsgs(tree, :x3)
dwnmsg = getDwnMsgs(tree, :x1)
dwnmsg = getDwnMsgs(tree, :x0)

printCliqHistorySummary(tree, :x1)

## test second inference/init clique for correct inferred dimensions


printCliqHistorySummary(tree, :x2)
getCliq(tree, :x2)

sfg = hist[5][17][4].cliqSubFg

isInitialized(sfg,:x1)
getVariableInferredDim(sfg,:x1)
isInitialized(fg,:x1)
getVariableInferredDim(fg,:x1)

isInitialized(sfg,:x3)
getVariableInferredDim(sfg,:x3)
isInitialized(fg,:x3)
getVariableInferredDim(fg,:x3)

isInitialized(sfg,:x2)
getVariableInferredDim(sfg,:x2)
isInitialized(fg,:x2)
getVariableInferredDim(fg,:x2)




sfg = hist[5][17][4].cliqSubFg
upmsg = prepCliqInitMsgsUp(sfg, getCliq(tree, :x2))
upMsg(tree, :x2)






## test root inference/init clique for correct inferred dimensions


printCliqHistorySummary(tree, :x3)
getCliq(tree, :x3)

sfg = hist[1][10][4].cliqSubFg

isInitialized(sfg,:x3)
getVariableInferredDim(sfg,:x3)
isInitialized(fg,:x3)
getVariableInferredDim(fg,:x3)

isInitialized(sfg,:x5)
getVariableInferredDim(sfg,:x5)
isInitialized(fg,:x5)
getVariableInferredDim(fg,:x5)

isInitialized(sfg,:l1)
getVariableInferredDim(sfg,:l1)
isInitialized(fg,:l1)
getVariableInferredDim(fg,:l1)


sfg = hist[1][10][4].cliqSubFg
upmsg = prepCliqInitMsgsUp(sfg, getCliq(tree, :x3))
upMsg(tree, :x3)







## test root inference/init clique for correct inferred dimensions


printCliqHistorySummary(tree, :x1)
getCliq(tree, :x1)

sfg = hist[2][16][4].cliqSubFg
getVariableInferredDim(sfg,:x1)

sfg = hist[2][17][4].cliqSubFg
getVariableInferredDim(sfg,:x1)



stuff = sandboxCliqResolveStep(tree,:x1,16)



isInitialized(fg,:x1)
getVariableInferredDim(fg,:x1)

isInitialized(sfg,:l1)
getVariableInferredDim(sfg,:l1)
isInitialized(fg,:l1)
getVariableInferredDim(fg,:l1)

isInitialized(sfg,:x3)
getVariableInferredDim(sfg,:x3)
isInitialized(fg,:x3)
getVariableInferredDim(fg,:x3)





sfg = hist[2][16][4].cliqSubFg
upmsg = prepCliqInitMsgsUp(sfg, getCliq(tree, :x1))
getUpMsgs(tree, :x3)








# make sure all infer dims are properly set

isInitialized(fg, :x0)
getVariableInferredDim(fg, :x0)

isInitialized(fg, :x1)
getVariableInferredDim(fg, :x1)

isInitialized(fg, :x2)
getVariableInferredDim(fg, :x2)

isInitialized(fg, :x3)
getVariableInferredDim(fg, :x3)

isInitialized(fg, :x4)
getVariableInferredDim(fg, :x4)

isInitialized(fg, :x5)
getVariableInferredDim(fg, :x5)

isInitialized(fg, :x6)
getVariableInferredDim(fg, :x6)

isInitialized(fg, :l1)
getVariableInferredDim(fg, :l1)






## Hex 2



posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, labels=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br, autoinit=false )


posecount = driveHex(fg, posecount)


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )

# writeGraphPdf(fg,show=true)


getSolverParams(fg).downsolve = false


tree, smt, hist = solveTree!(fg, tree, recordcliqs=ls(fg))










## debugging


# fetchCliqTaskHistoryAll!(smt, hist)
# assignTreeHistory!(tree, hist)
# printCliqHistorySummary(hist, tree, :x10)
# makeCsmMovie(fg, tree)



# down msg to :x10
dwinmsgs = prepCliqInitMsgsDown!(fg, tree, getCliq(tree, :x9), getCliq(tree,:x11), dbgnew=true)



# up msg to :x6
getCliq(tree, :x6)
printCliqHistorySummary(hist, tree, :x6)
sfg = hist[11][9][4].cliqSubFg
# writeGraphPdf(sfg, show=true)
upmsg = prepCliqInitMsgsUp(sfg, getCliq(tree, :x6))


## PROBLEM: initializes and solves, but inferdim not populated
# deeper issue is that neither x5 or x6 have valid inferdim values
sfg = hist[11][8][4].cliqSubFg
isInitialized(sfg,:x7)
getVariableInferredDim(sfg,:x7)
isInitialized(sfg,:x5)
getVariableInferredDim(sfg,:x5)
sfg = hist[11][9][4].cliqSubFg
isInitialized(sfg,:x7)
getVariableInferredDim(sfg,:x7)



cliq = hist[11][8][4].cliq
put!(getData(cliq).initUpChannel, :upsolved)

take!(getData(cliq).initUpChannel)

stuff = sandboxCliqResolveStep(tree,:x6,8)









#
