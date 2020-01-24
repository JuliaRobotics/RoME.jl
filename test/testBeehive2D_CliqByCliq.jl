
# using Revise

using Test

using RoME


include(joinpath(@__DIR__,"BeehiveTestUtils.jl"))


@testset "sanity check on Hex example" begin

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

getSolverParams(fg).drawtree = false
getSolverParams(fg).showtree = false
getSolverParams(fg).downsolve = false
getSolverParams(fg).multiproc = false
getSolverParams(fg).async = false
# to disable parent factor sharing
# getSolverParams(fg).devParams[:dontUseParentFactorsInitDown] = ""



# direct solve would be
tree, smt, hist = solveTree!(fg) #, recordcliqs=ls(fg))


@test 45 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[1,:] .< 3.0)
@test 45 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[2,:] .< 3.0)
@test 45 < sum(-0.3 .< getPoints(getKDE(fg, :x0))[3,:] .< 0.3)

@test 45 < sum(7.0 .< getPoints(getKDE(fg, :x1))[1,:] .< 13.0)
@test 45 < sum(-3.0 .< getPoints(getKDE(fg, :x1))[2,:] .< 3.0)
@test 45 < sum(0.7 .< getPoints(getKDE(fg, :x1))[3,:] .< 1.3)

@test 45 < sum(12.0 .< getPoints(getKDE(fg, :x2))[1,:] .< 18.0)
@test 45 < sum(6.0 .< getPoints(getKDE(fg, :x2))[2,:] .< 11.0)
@test 45 < sum(1.8 .< getPoints(getKDE(fg, :x2))[3,:] .< 2.4)

@test 45 < sum(7.0 .< getPoints(getKDE(fg, :x3))[1,:] .< 13.0)
@test 45 < sum(15.0 .< getPoints(getKDE(fg, :x3))[2,:] .< 20.0)
# @test 45 < sum(-0.3 .< getPoints(getKDE(fg, :x3))[3,:] .< 0.3)

@test 45 < sum(-5.0 .< getPoints(getKDE(fg, :x4))[1,:] .< 5.0)
@test 45 < sum(13.0 .< getPoints(getKDE(fg, :x4))[2,:] .< 22.0)
@test 45 < sum(-2.8 .< getPoints(getKDE(fg, :x4))[3,:] .< -1.5)

@test 45 < sum(-8.0 .< getPoints(getKDE(fg, :x5))[1,:] .< -2.0)
@test 45 < sum(6.0 .< getPoints(getKDE(fg, :x5))[2,:] .< 11.0)
@test 45 < sum(-1.3 .< getPoints(getKDE(fg, :x5))[3,:] .< -0.7)

@test 45 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[1,:] .< 3.0)
@test 45 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[2,:] .< 3.0)
@test 45 < sum(-0.3 .< getPoints(getKDE(fg, :x6))[3,:] .< 0.3)

@test 45 < sum(17.0 .< getPoints(getKDE(fg, :l1))[1,:] .< 23.0)
@test 45 < sum(-5.0 .< getPoints(getKDE(fg, :l1))[2,:] .< 5.0)



# #  Do some plotting
# using RoMEPlotting
#
# Gadfly.set_default_plot_size(35cm,20cm)
# pl = drawPosesLandms(fg, meanmax=:max, drawhist=true)
# pl |> PDF("/tmp/caesar/test.pdf");  @async run(`evince /tmp/caesar/test.pdf`)
# pl |> PNG("/tmp/caesar/test.png");
# drawTree(tree, imgs=true)


end # testset



## OR manually do clique by clique
#
# # solve by hand, one cliq at a time
# tree = wipeBuildNewTree!(fg)
# # drawTree(tree,show=true)
#
#
# # solve the first cliq
# smt, hist = solveCliq!(fg, tree, :x0, recordcliq=true)
#
# # solve second clq
# smt, hist = solveCliq!(fg, tree, :x2, cliqHistories=hist, recordcliq=true)
# # resetCliqSolve!(fg, tree, :x2)
#
# # solve forth clq
# smt, hist = solveCliq!(fg, tree, :x4, cliqHistories=hist, recordcliq=true)
# # solve forth clq
# smt, hist = solveCliq!(fg, tree, :x6, cliqHistories=hist, recordcliq=true)
# # solve forth clq
# smt, hist = solveCliq!(fg, tree, :x3, cliqHistories=hist, recordcliq=true)
