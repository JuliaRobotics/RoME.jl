
using Revise

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


function offsetHexLeg(dfg::G, posecount::Int; direction=:right) where G <: AbstractDFG
    psym = Symbol("x$(posecount-1)")
    nsym = Symbol("x$(posecount)")
    posecount += 1
    addVariable!(dfg, nsym, Pose2)
    pp = nothing
    if direction == :right
        pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    elseif direction == :left
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    end
    addFactor!(dfg, [psym; nsym], pp, autoinit=false )
    return posecount
end



## start with an empty factor graph object
fg = initfg()

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


# writeGraphPdf(fg,engine="neato")

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
getSolverParams(fg).downsolve = false
getSolverParams(fg).multiproc = false
getSolverParams(fg).async = true
getSolverParams(fg).devParams[:useParentFactorsInitDown] = ""



# direct solve would be
# tree, smt, hist = solveTree!(fg)


# tree, smt, hist = solveTree!(fg, skipcliqids=[:x1;:x6;:x4;:x3], recordcliqs=[:x2;:x0])


# solve by hand, one cliq at a time
tree = wipeBuildNewTree!(fg)
# drawTree(tree,show=true)


# solve the first cliq
smt, hist = solveCliq!(fg, tree, :x0, recordcliq=true)


# solve second clq
smt, hist = solveCliq!(fg, tree, :x2, cliqHistories=hist, recordcliq=true)
# resetCliqSolve!(fg, tree, :x2)



# DEBUG after solveCliq :x2


# solve third clq
smt, hist = solveCliq!(fg, tree, :x1, cliqHistories=hist, recordcliq=true)


# Plot to see what is going on
Gadfly.set_default_plot_size(35cm,20cm)
plotKDE(fg, [:x0;:x1;:x2], dims=[1;2],levels=1)

plotKDE(fg, :x1, dims=[1;2],levels=1)





## DEBUG where is down message :x3

prnt = getCliq(tree, :x3)
getCliqInitUpMsgs(prnt)


assignTreeHistory!(tree, hist)
printCliqHistorySummary(tree, :x1)


# think the issue is in here
# 21:01:00.856  8   null         attemptCliqInitUp     false null | upsolved upsolved



# somethings up with cliq x1

csmc_8_test = getCliqSolveHistory(tree, :x1)[8][4]
csmc_9_test = getCliqSolveHistory(tree, :x1)[9][4]

# Just before the mistake
drawTree(csmc_8_test.tree, show=true)
prnt_8_test = getCliq(csmc_8_test.tree, :x3)
getCliqInitUpMsgs(prnt_8_test)

# just after the mistake
prnt_9_test = getCliq(csmc_9_test.tree, :x3)
getCliqInitUpMsgs(prnt_9_test)


# so develop  based on sandbox step 8
stuff = sandboxCliqResolveStep(tree,:x1,8)
getCliqInitUpMsgs(getCliq(stuff[4].tree,:x3))






## DEBUG





# solve forth clq
smt, hist = solveCliq!(fg, tree, :x4, cliqHistories=hist, recordcliq=true)
# solve forth clq
smt, hist = solveCliq!(fg, tree, :x6, cliqHistories=hist, recordcliq=true)
# solve forth clq
smt, hist = solveCliq!(fg, tree, :x3, cliqHistories=hist, recordcliq=true)




# tree, smt, chi = solveTree!(fg, recordcliqs=ls(fg))

# hist = getCliqSolveHistory(tree, :x1)


Gadfly.set_default_plot_size(35cm,20cm)
drawPosesLandms(fg, meanmax=:max) |> PDF("/tmp/test.pdf");  @async run(`evince /tmp/test.pdf`)
# tree = wipeBuildNewTree!(fg)
# drawTree(tree, imgs=true)



#
#
# plotTreeUpMsgs(fg, tree, :x1, levels=1)
# plotTreeUpMsgs(fg, tree, :x3, levels=1)
# plotTreeUpMsgs(fg, tree, :x5, levels=1)
# plotTreeUpMsgs(fg, tree, :l1, levels=1)






cliq = getCliq(tree, :x2)

# OLD
dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
plotKDE(dwinmsgs[:x1][1], dims=[1;2], levels=2)

# NEW
dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
plotKDE(dwinmsgs[:x1][1], dims=[1;2], levels=2)








## Check init message for x3

cliq = getCliq(tree, :x1)

# OLD
dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)

# NEW
dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)





## Check init message for x4

cliq = getCliq(tree, :x4)

# OLD
dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)

# NEW
dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)







plotKDE(fg, :x1, dims=[1;2])
plotPose(fg, :x1)

plotKDE(fg, [:x0;:x1;:x2;:x4;:x6], dims=[1;2],levels=1)


drawTree(tree, imgs=true)




## check contents

cliq = getCliq(tree, :x1)

getCliqMsgsUp(cliq)


#











## DEBUG after solveCliq :x2 difference in down init cycle when :useParentFactorsInitDown

assignTreeHistory!(tree, hist)
printCliqHistorySummary(tree, :x2)

# WORKS

delete!(getSolverParams(getCliqSolveHistory(tree, :x2)[15][4].dfg).devParams, :useParentFactorsInitDown)

good = sandboxCliqResolveStep(tree,:x2,15)

# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x1, title="(X2,) sandbox 15 good,")
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x2, title="(X2,) sandbox 15 good,")


# BREAKS

getSolverParams(getCliqSolveHistory(tree, :x2)[15][4].dfg).devParams[:useParentFactorsInitDown] = ""

bad = sandboxCliqResolveStep(tree,:x2,15)

# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x1, title="(X2,) sandbox 15 bad,")
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x2, title="(X2,) sandbox 15 bad,")

plotPairVariables(good[4].cliqSubFg, bad[4].cliqSubFg, :x2, title="(X2,) sandbox 15,")
# plotPairPose2(good[4].cliqSubFg, bad[4].cliqSubFg, :x2, title="(X2,) sandbox 15,")

Gadfly.set_default_plot_size(30cm,18cm)

## DEBUG













##



#3 First difference, up msg from (:x2) cliq is 'wrong'

Gadfly.set_default_plot_size(30cm,18cm)
um1 = IIF.getCliqMsgsUp(tree, :x0)
plotKDE(um1[:x1], dims=[1;2], levels=2)



um2 = IIF.getCliqMsgsUp(tree, :x2)
plotKDE(um2[:x1], dims=[1;2], levels=2)




csmcX2_test = getCliqSolveHistory(tree, :x2)

# csmcX2_test[15][4].cliqSubFg

plotKDE(getKDE(csmcX2_test[17][4].cliqSubFg, :x1),dims=[1;2],levels=2,title="csmcX2_test[17]")






## SOMETHING changes between csmcX2_ref[16/17] and csmcX2_test[16/17]

writeGraphPdf(csmcX2_ref[17][4].cliqSubFg, show=true)
writeGraphPdf(csmcX2_test[17][4].cliqSubFg, show=true)


plotKDE(csmcX2_ref[16][4].cliqSubFg, [:x1;:x2;:x3], dims=[1;2])
plotKDE(csmcX2_test[16][4].cliqSubFg, [:x1;:x2;:x3], dims=[1;2])



isInitialized(csmcX2_test[16][4].cliqSubFg, :x3)








# before init cycle
plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x1)
plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x2)
plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x3)


# after solve
plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x1)
plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x2)
plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x3)
