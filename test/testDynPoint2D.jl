# test the DynPoint2 functionality

using RoME, IncrementalInference, Distributions
using Base: Test


# N = 75
fg = initfg()

# add two point locations
v0 = addNode!(fg, :x0, DynPoint2(ut=0))

v1 = addNode!(fg, :x1, DynPoint2(ut=1000_000))

# Prior factor as boundary condition
pp0 = DynPoint2VelocityPrior(MvNormal([zeros(2);10*ones(2)], 0.1*eye(4)))
f0 = addFactor!(fg, [:x0;], pp0)

# conditional likelihood between Dynamic Point2
dp2dp2 = DynPoint2DynPoint2(MvNormal([10*ones(2);zeros(2)], 0.1*eye(4)))
f1 = addFactor!(fg, [:x0;:x1], dp2dp2)

# Graphs.plot(fg.g)
ensureAllInitialized!(fg)




tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



# using RoMEPlotting,
using KernelDensityEstimate #, KernelDensityEstimatePlotting

# X1 = getVal(fg, :x1)
@show x0 = getKDEMax(getVertKDE(fg, :x0))

@show x1 = getKDEMax(getVertKDE(fg, :x1))





# N = 75
fg = initfg()

# add two point locations
v0 = addNode!(fg, :x0, DynPoint2(ut=0))
v1 = addNode!(fg, :x1, DynPoint2(ut=1000_000))
v2 = addNode!(fg, :x2, DynPoint2(ut=2000_000))


# Prior factor as boundary condition
pp0 = DynPoint2VelocityPrior(MvNormal([zeros(2);10*ones(2)], 0.1*eye(4)))
f0 = addFactor!(fg, [:x0;], pp0)

# conditional likelihood between Dynamic Point2
dp2dp2 = VelPoint2VelPoint2(MvNormal([10*ones(2);zeros(2)], 0.1*eye(4)))
f1 = addFactor!(fg, [:x0;:x1], dp2dp2)

# conditional likelihood between Dynamic Point2
dp2dp2 = VelPoint2VelPoint2(MvNormal([10*ones(2);zeros(2)], 0.1*eye(4)))
f2 = addFactor!(fg, [:x1;:x2], dp2dp2)


# Graphs.plot(fg.g)
ensureAllInitialized!(fg)




tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



# X1 = getVal(fg, :x1)
@show x0 = getKDEMax(getVertKDE(fg, :x0))

@show x1 = getKDEMax(getVertKDE(fg, :x1))

@show x2 = getKDEMax(getVertKDE(fg, :x2))











@testset "test VelPoint2VelPoint2" begin


N=100
pμ = [0.0,0,1,0]
pσ = diagm([0.1;0.1;0.1;0.1].^2)

fg = initfg();

addNode!(fg, :x1, DynPoint2(ut=0))
pp = DynPoint2VelocityPrior(MvNormal(pμ,pσ))
addFactor!(fg, [:x1;], pp, autoinit=false)

addNode!(fg, :x2, DynPoint2(ut=1_000_000))
dpμ = [1.0;0;0;0];
dpσ = diagm([1.0;1;0.5;0.01].^2)
vp = VelPoint2VelPoint2(MvNormal(dpμ,dpσ))
addFactor!(fg, [:x1,:x2], vp, autoinit=false)

# writeGraphPdf(fg)

# note, nothign is initialized yet...
isInitialized(fg, :x1)
isInitialized(fg, :x2)

# lets init the first
doautoinit!(fg, :x1)

X1 = getVal(fg, :x1)
@test 0.7*N < sum(-0.5 .< X1[1,:] .< 0.5)
@test 0.7*N < sum(-0.5 .< X1[2,:] .< 0.5)

# now test "projection" / probabilistic-convolution through VelPointVelPoint{D, N}
pX2 = approxConv(fg, :x1x2f1, :x2)

@test 0.7*N < sum(0 .< pX2[1,:] .< 2)
@test 0.7*N < sum(-1 .< pX2[2,:] .< 1)


# using RoMEPlotting, Gadfly
# using KernelDensityEstimate
# using KernelDensityEstimatePlotting
# X1 = getVal(fg,:x1)
# plot(x=pX2[1,:], y=pX2[2,:], Geom.hexbin)
# pl = plotKDE(kde!(pX2[1:2,:]), dims=[1;2], levels=3, c=["blue"])
# pl = plotKDE(kde!(X1[1:2,:]), dims=[1;2], levels=3, c=["blue"])

end






#
