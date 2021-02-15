# Taken from #380
# added during development for JuliaRobotics/IncrementalInference.jl#1051


using Test
using RoME
using Statistics


##

@testset "test bearing range with inflation, #380, IIF #1051" begin

##

fg = initfg()

getSolverParams(fg).inflation = 3.0

pr_noise = [0.01, 0.01, 0.001]
od_noise = [0.2;0.2;0.2]
σ_bearing = 0.008
σ_range = 0.01

addVariable!(fg, :x0, Pose2)
prpo = MvNormal([0,0,-pi], pr_noise)
addFactor!(fg, [:x0], PriorPose2(MvNormal(rand(prpo), pr_noise)))

addVariable!(fg, :l1, Point2, tags=[:LANDMARK])
addVariable!(fg, :l2, Point2, tags=[:LANDMARK])
p2br = Pose2Point2BearingRange(Normal(pi/4 + rand(Normal(0,σ_bearing)), σ_bearing),
                                Normal(sqrt(2) + rand(Normal(0,σ_range)), σ_range))
addFactor!(fg, [:x0; :l1], p2br)

addVariable!(fg, :x1, Pose2)
pp = MvNormal([1.0,0,0], od_noise)
addFactor!(fg, [:x0,:x1], Pose2Pose2(MvNormal(rand(pp), od_noise)))

p2br = Pose2Point2BearingRange(Normal(pi/2 + rand(Normal(0,σ_bearing)), σ_bearing),
                                Normal(1 + rand(Normal(0,σ_range)), σ_range))
addFactor!(fg, [:x1; :l1], p2br)

p2br = Pose2Point2BearingRange(Normal(-pi/2 + rand(Normal(0,σ_bearing)), σ_bearing),
                                Normal(1 + rand(Normal(0,σ_range)), σ_range))
addFactor!(fg, [:x1; :l2], p2br)

addVariable!(fg, :x2, Pose2)
pp = MvNormal([1.0,0,0], od_noise)
addFactor!(fg, [:x1,:x2], Pose2Pose2(MvNormal(rand(pp), od_noise)))

p2br = Pose2Point2BearingRange(Normal(3pi/4 + rand(Normal(0,σ_bearing)), σ_bearing),
                                Normal(sqrt(2) + rand(Normal(0,σ_range)), σ_range))
addFactor!(fg, [:x2; :l1], p2br)

p2br = Pose2Point2BearingRange(Normal(-3pi/4 + rand(Normal(0,σ_bearing)), σ_bearing),
                                Normal(sqrt(2) + rand(Normal(0,σ_range)), σ_range))
addFactor!(fg, [:x2; :l2], p2br)


##

# nonparametric solution
solveGraph!(fg);
IIF.solveFactorGraphParametric!(fg)


##

@show getPPE(fg, :x2, :parametric).suggested
@show getPPE(fg, :x2, :default).suggested

test_err = getPPE(fg, :x2, :default).suggested - getPPE(fg, :x2, :parametric).suggested
# arg, workaround until #244
theta = (getBelief(fg, :x2, :default) |> getPoints)[3,:] .+ pi
theta .= TU.wrapRad.(theta)
theta_ = Statistics.mean(theta) + pi
test_err[3] = abs(getPPE(fg, :x2, :parametric).suggested[3]) - theta_
@show test_err .= abs.(test_err)

@test isapprox(test_err[1], 0, atol=0.5)
@test isapprox(test_err[2], 0, atol=0.5)
@test isapprox(test_err[3], 0, atol=0.5)


##

end


## 

# using RoMEPlotting
# Gadfly.set_default_plot_size(25cm,20cm)

# ##

# pl1 = plotSLAM2D(fg, solveKey=:default, drawPoints=true, drawEllipse=true, drawContour=false, xmin=-2.5,xmax=0,ymin=-1.5,ymax=1.5)
# pl2 = plotSLAM2D(fg, solveKey=:parametric, drawPoints=false, drawContour=false, xmin=-2.5,xmax=0,ymin=-1.5,ymax=1.5)
# vstack(pl1, pl2)

##

# plotPose(fg, :x2)

##