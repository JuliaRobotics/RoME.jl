
# using Revise
# using Distributed
# addprocs(4)

using RoME
# @everywhere using RoME

using Test


include(joinpath(@__DIR__,"BeehiveTestUtils.jl"))


@testset "sanity check on Hex example" begin

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
getSolverParams(fg).multiproc = false
# fg.solverParams.async = true

# solve
tree, smt, chi = solveTree!(fg, recordcliqs=ls(fg))




## See visualization

# hist = getCliqSolveHistory(tree, :x1)

#  Do some plotting
using RoMEPlotting
Gadfly.set_default_plot_size(35cm,25cm)
drawPosesLandms(fg, meanmax=:max) |> PDF("/tmp/test.pdf");  @async run(`evince /tmp/test.pdf`)




## check outcome on first hex

@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[1,:] .< 3.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[2,:] .< 3.0)
@test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x0))[3,:] .< 0.3)

@test 80 < sum(7.0 .< getPoints(getKDE(fg, :x1))[1,:] .< 13.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x1))[2,:] .< 3.0)
@test 80 < sum(0.7 .< getPoints(getKDE(fg, :x1))[3,:] .< 1.3)

@test 80 < sum(12.0 .< getPoints(getKDE(fg, :x2))[1,:] .< 18.0)
@test 80 < sum(6.0 .< getPoints(getKDE(fg, :x2))[2,:] .< 11.0)
@test 80 < sum(1.8 .< getPoints(getKDE(fg, :x2))[3,:] .< 2.4)

@test 80 < sum(7.0 .< getPoints(getKDE(fg, :x3))[1,:] .< 13.0)
@test 80 < sum(15.0 .< getPoints(getKDE(fg, :x3))[2,:] .< 20.0)
# @test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x3))[3,:] .< 0.3)

@test 80 < sum(-4.0 .< getPoints(getKDE(fg, :x4))[1,:] .< 4.0)
@test 80 < sum(15.0 .< getPoints(getKDE(fg, :x4))[2,:] .< 20.0)
@test 80 < sum(-2.4 .< getPoints(getKDE(fg, :x4))[3,:] .< -1.8)

@test 80 < sum(-8.0 .< getPoints(getKDE(fg, :x5))[1,:] .< -2.0)
@test 80 < sum(6.0 .< getPoints(getKDE(fg, :x5))[2,:] .< 11.0)
@test 80 < sum(-1.3 .< getPoints(getKDE(fg, :x5))[3,:] .< -0.7)

@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[1,:] .< 3.0)
@test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[2,:] .< 3.0)
@test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x6))[3,:] .< 0.3)

@test 80 < sum(17.0 .< getPoints(getKDE(fg, :l1))[1,:] .< 23.0)
@test 80 < sum(-5.0 .< getPoints(getKDE(fg, :l1))[2,:] .< 5.0)



# check outcome on second hex



# @test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[1,:] .< 3.0)
# @test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x0))[2,:] .< 3.0)
# @test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x0))[3,:] .< 0.3)
#
# @test 80 < sum(7.0 .< getPoints(getKDE(fg, :x1))[1,:] .< 13.0)
# @test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x1))[2,:] .< 3.0)
# @test 80 < sum(0.7 .< getPoints(getKDE(fg, :x1))[3,:] .< 1.3)
#
# @test 80 < sum(12.0 .< getPoints(getKDE(fg, :x2))[1,:] .< 18.0)
# @test 80 < sum(6.0 .< getPoints(getKDE(fg, :x2))[2,:] .< 11.0)
# @test 80 < sum(1.8 .< getPoints(getKDE(fg, :x2))[3,:] .< 2.4)
#
# @test 80 < sum(7.0 .< getPoints(getKDE(fg, :x3))[1,:] .< 13.0)
# @test 80 < sum(15.0 .< getPoints(getKDE(fg, :x3))[2,:] .< 20.0)
# # @test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x3))[3,:] .< 0.3)
#
# @test 80 < sum(-4.0 .< getPoints(getKDE(fg, :x4))[1,:] .< 4.0)
# @test 80 < sum(15.0 .< getPoints(getKDE(fg, :x4))[2,:] .< 20.0)
# @test 80 < sum(-2.4 .< getPoints(getKDE(fg, :x4))[3,:] .< -1.8)
#
# @test 80 < sum(-8.0 .< getPoints(getKDE(fg, :x5))[1,:] .< -2.0)
# @test 80 < sum(6.0 .< getPoints(getKDE(fg, :x5))[2,:] .< 11.0)
# @test 80 < sum(-1.3 .< getPoints(getKDE(fg, :x5))[3,:] .< -0.7)
#
# @test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[1,:] .< 3.0)
# @test 80 < sum(-3.0 .< getPoints(getKDE(fg, :x6))[2,:] .< 3.0)
# @test 80 < sum(-0.3 .< getPoints(getKDE(fg, :x6))[3,:] .< 0.3)
#
# @test 80 < sum(17.0 .< getPoints(getKDE(fg, :l1))[1,:] .< 23.0)
# @test 80 < sum(-5.0 .< getPoints(getKDE(fg, :l1))[2,:] .< 5.0)




end # testset
