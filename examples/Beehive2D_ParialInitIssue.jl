
# using Revise
# using Distributed
# addprocs(4)

using RoME
# @everywhere using RoME

#  Do some plotting
using RoMEPlotting
# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

0

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

# fg.solverParams
# fg.solverParams.isfixedlag = true
# fg.solverParams.qfl = 20
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



## hex 2

posecount = offsetHexLeg(fg, posecount, direction=:right)

# Add landmarks with Bearing range measurements
addVariable!(fg, :l2, Point2, labels=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br, autoinit=false )


posecount = driveHex(fg, posecount, steps=5)



# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [Symbol("x$(posecount-1)"); :l2], p2br2, autoinit=false )



# writeGraphPdf(fg,engine="neato")

getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
getSolverParams(fg).downsolve = false
fg.solverParams.async = true


tree, smt, chi = solveTree!(fg, recordcliqs=ls(fg))

# hist = getCliqSolveHistory(tree, :x1)



Gadfly.set_default_plot_size(35cm,25cm)
drawPosesLandms(fg, meanmax=:max) |> PDF("/tmp/test.pdf");  @async run(`evince /tmp/test.pdf`)
# tree = wipeBuildNewTree!(fg)
# drawTree(tree, imgs=true)

#
# using Dates
#
# dd = [true;]
# mkpath("/tmp/caesar/anitree")
# @async begin
#   global dd
#   i =0
#   while dd[1]
#     i += 1
#     tt = split(string(now()), 'T')[end]
#     drawTree(tree, filepath="/tmp/caesar/anitree/bt_$i.pdf", show=false)
#     sleep(0.1)
#   end
# end
#
# # dd[1] = false


#
