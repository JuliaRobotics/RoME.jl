using RoME
using Statistics
# , Distributions
using Test

@testset "test sampling from BearingRange factor..." begin

global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))

global meas = getSample(p2br, 100)
@test abs(Statistics.mean(meas[1][1,:])) < 0.1
@test 0.05 < abs(Statistics.std(meas[1][1,:])) < 0.2

@test abs(Statistics.mean(meas[1][2,:]) - 20.0) < 1.0
@test 0.5 < abs(Statistics.std(meas[1][2,:])) < 1.5

end


@testset "test BearingRange factor residual function..." begin

global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))

global xi = zeros(3,1)
global li = zeros(2,1); li[1,1] = 20.0;
global zi = (zeros(2,1),); zi[1][2,1] = 20.0

global idx = 1
global res = zeros(2)
p2br(res, IncrementalInference.FactorMetadata(), idx, zi, xi, li)
@show res

@test norm(res) < 1e-14


global xi = zeros(3,1)
global li = zeros(2,1); li[2,1] = 20.0;
global zi = (zeros(2,1),); zi[1][:,1] = [pi/2;20.0]

global idx = 1
global res = zeros(2)
p2br(res, IncrementalInference.FactorMetadata(), idx, zi, xi, li)
@show res
@test norm(res) < 1e-14


global xi = zeros(3,1); xi[3,1] = pi/2
global li = zeros(2,1); li[2,1] = 20.0;
global zi = (zeros(2,1),); zi[1][:,1] = [0.0;20.0]

global idx = 1
global res = zeros(2)
p2br(res, IncrementalInference.FactorMetadata(), idx, zi, xi, li)
@show res
@test norm(res) < 1e-14


global xi = zeros(3,2); xi[3,2] = -pi/2
global li = zeros(2,2); li[1,2] = 20.0;
# zi = ([0.0;pi/2],[0.0;20.0],)
global zi = (zeros(2,2),); zi[1][:,2] = [pi/2;20.0]

global idx = 2
global res = zeros(2)
p2br(res, IncrementalInference.FactorMetadata(), idx, zi, xi, li)
@show res
@test norm(res) < 1e-14

end




@testset "test unimodal bearing range factor, solve for landmark..." begin

# Start with an empty graph
global N = 1
global fg = initfg()

#add pose with partial constraint
addNode!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I, 3,3))), autoinit=false) # TODO IIF.Prior
# force particular initialization
setVal!(fg, :x0, zeros(3,1))

##----------- sanity check that predictbelief plumbing is doing the right thing
global pts = predictbelief(fg, :x0, ls(fg, :x0), N=75)
@test sum(abs.(Statistics.mean(pts,dims=2)) .< [0.1; 0.1; 0.1]) == 3
@test sum([0.05; 0.05; 0.05] .< Statistics.std(pts,dims=2) .< [0.15; 0.15; 0.15]) == 3
#------------

# Add landmark
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
global li = zeros(2,1); li[1,1] = 20.0;
setVal!(fg, :l1, li)


# Add bearing range measurement between pose and landmark
global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false)

# there should be just one (the bearingrange) factor connected to :l1
@test length(ls(fg, :l1)) == 1
# writeGraphPdf(fg)

# check the forward convolution is working properly
global pts = predictbelief(fg, :l1, ls(fg, :l1), N=75)
@test sum(abs.(Statistics.mean(pts,dims=2) - [20.0; 0.0]) .< [2.0; 2.0]) == 2
@test sum([0.1; 0.1] .< Statistics.std(pts,dims=2) .< [3.0; 3.0]) == 2

# using Gadfly, KernelDensityEstimate, KernelDensityEstimatePlotting
#
# pl = plotKDE(kde!(pts))
# pl.coord = Coord.Cartesian(xmin=-5,xmax=25, ymin=-10.0,ymax=10)
# pl

end


@testset "test unimodal bearing range factor, solve for pose..." begin

# Start with an empty graph
global N = 1
global fg = initfg()

# Add landmark
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
addFactor!(fg, [:l1], PriorPoint2(MvNormal([20.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))),  autoinit=false ) # TODO IIF.Prior
global li = zeros(2,1); li[1,1] = 20.0;
setVal!(fg, :l1, li)

#add pose with partial constraint
addNode!(fg, :x0, Pose2)
# force particular initialization
setVal!(fg, :x0, zeros(3,1))

# Add bearing range measurement between pose and landmark
global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false)

# there should be just one (the bearingrange) factor connected to :l1
@test length(ls(fg, :x0)) == 1
# writeGraphPdf(fg)

# check the forward convolution is working properly
global pts = predictbelief(fg, :x0, ls(fg, :x0), N=75)

@show abs.(Statistics.mean(pts,dims=2))
@test sum(abs.(Statistics.mean(pts,dims=2)) .< [2.0; 2.0; 2.0]) == 3
@show Statistics.std(pts,dims=2)
@test sum([0.1; 0.1; 0.01] .< Statistics.std(pts,dims=2) .< [5.0; 5.0; 2.0]) == 3

end
