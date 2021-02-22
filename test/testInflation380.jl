# Taken from #380
# added during development for JuliaRobotics/IncrementalInference.jl#1051


using Test
using RoME
using Statistics


##


@testset "test inflation is working via distance" begin

##

fg = initfg()
# getSolverParams(fg).inflation = 50.0

N = 100

addVariable!(fg, :x0, ContinuousEuclid{2})
addVariable!(fg, :x1, ContinuousEuclid{2})

X0_ = zeros(2,N)
X0_[1,:] .+= 1000.0
initManual!(fg, :x0, X0_)

addFactor!(fg, [:x0;:x1], EuclidDistance(Normal(100.0, 1.0)))

pts = approxConv(fg, :x0, :x1)


##
# does this give a "donut ring" at 1000?

res = 99999*ones(100)

for i in 1:N
  res[i] = calcFactorResidual(fg, :x0x1f1, [100.0], [1000.0;0.0], pts[:,i])[1]
end

@test 0.9*N < sum(abs.(res) .< 5)

## new test trying to force inflation error

X1_ = randn(2,N)
X1_[1,:] .+= 1100.0
initManual!(fg, :x1, X1_)


##

IIF._getCCW(fg, :x0x1f1).inflation = 50.0
pts = approxConv(fg, :x0x1f1, :x1)

initManual!(fg, :x1, pts)

pts = approxConv(fg, :x0x1f1, :x1)


##
# does this give a "donut ring" at 1000?

res = 99999*ones(100)

for i in 1:N
  res[i] = calcFactorResidual(fg, :x0x1f1, [100.0], [1000.0;0.0], pts[:,i])[1]
end

@test 0.9*N < sum(abs.(res) .< 5)


##

# using RoMEPlotting
# Gadfly.set_default_plot_size(25cm,20cm)

# ##

# # pts = getBelief(fg, :x1) |> getPoints
# plotKDE(manikde!(pts, ContinuousEuclid{2}))


##

end


@testset "test inflation on range solve" begin

##

N = 100
fg = initfg()
fg.solverParams.inflation = 10.0 # super size the inflation to force wide coverage

addVariable!(fg, :x1, ContinuousEuclid{2})
addVariable!(fg, :l1, ContinuousEuclid{2})

addFactor!(fg, [:l1], Prior(MvNormal([-1000.0,0], [0.1, 0.1])))
addFactor!(fg, [:x1; :l1], EuclidDistance(Normal(100.0, 1.0)))

pts = zeros(2,100)
pts[1,:] .-= 900
initManual!(fg, :x1, pts)

##

tree, _, = solveGraph!(fg);

##

pts = getBelief(fg, :x1) |> getPoints

@test 0.9*N < sum( -1150 .< pts[1,:] .< -850)
@test 0.9*N < sum( -150 .< pts[2,:] .< 150)

pts_ = [norm(pts[:,i] - [-1000;0]) for i in 1:N]

@test 0.9*N < sum(80 .< pts_ .< 120)

# must still test spread

@test 0.2*N < sum(-1150 .< pts[1,:] .< -1000)
@test 0.2*N < sum(-1000 .< pts[1,:] .< -850)
@test 0.2*N < sum(-150 .< pts[2,:] .< 0)
@test 0.2*N < sum(0 .< pts[2,:] .< 150)


##

# pl = plotKDE(fg, ls(fg))

##

end
  

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
@show theta_ = Statistics.mean(theta)
@show theta_ref = TU.wrapRad(getPPE(fg, :x2, :parametric).suggested[3] + pi)
test_err[3] = theta_ref - theta_
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