
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


# direct solve would be
# tree, smt, hist = solveTree!(fg)

# solve by hand, one cliq at a time
tree = wipeBuildNewTree!(fg, drawpdf=true)


# solve the first cliq
smt, hist = solveCliq!(fg, tree, :x0)
# solve second clq
smt, hist = solveCliq!(fg, tree, :x2, cliqHistories=hist)



#3 First difference, up msg from (:x2) cliq is 'wrong'

# Gadfly.set_default_plot_size(35cm,20cm)
um1 = IIF.getCliqMsgsUp(tree, :x0)
plotKDE(um1[:x1], dims=[1;2], levels=2)



um2 = IIF.getCliqMsgsUp(tree, :x2)
plotKDE(um2[:x1], dims=[1;2], levels=2)



resetCliqSolve!(fg, tree, :x2)

drawTree(tree)







# solve third clq
smt, hist = solveCliq!(fg, tree, :x1, cliqHistories=hist)


# Plot to see what is going on
Gadfly.set_default_plot_size(35cm,25cm)
plotKDE(fg, [:x0;:x1;:x2], dims=[1;2],levels=1)



# solve forth clq
smt, hist = solveCliq!(fg, tree, :x4, cliqHistories=hist)
# solve forth clq
smt, hist = solveCliq!(fg, tree, :x6, cliqHistories=hist)
# solve forth clq
smt, hist = solveCliq!(fg, tree, :x3, cliqHistories=hist)




# tree, smt, chi = solveTree!(fg, recordcliqs=ls(fg))

# hist = getCliqSolveHistory(tree, :x1)


Gadfly.set_default_plot_size(35cm,25cm)
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
