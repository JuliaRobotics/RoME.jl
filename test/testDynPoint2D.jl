# test the DynPoint2 functionality

using RoME
using Test
# using RoMEPlotting,
# , IncrementalInference, Distributions
# using KernelDensityEstimate #, KernelDensityEstimatePlotting



@testset "test DynPoint2..." begin

# N = 75
global fg = initfg()

# add two point locations
global v0 = addNode!(fg, :x0, DynPoint2(ut=0))

global v1 = addNode!(fg, :x1, DynPoint2(ut=1000_000))

# Prior factor as boundary condition
global pp0 = DynPoint2VelocityPrior(MvNormal([zeros(2);10*ones(2)], 0.1*Matrix{Float64}(LinearAlgebra.I, 4,4)))
global f0 = addFactor!(fg, [:x0;], pp0)
# conditional likelihood between Dynamic Point2
global dp2dp2 = DynPoint2DynPoint2(MvNormal([10*ones(2);zeros(2)], 0.1*Matrix{Float64}(LinearAlgebra.I, 4,4)))
global f1 = addFactor!(fg, [:x0;:x1], dp2dp2)

# Graphs.plot(fg.g)
ensureAllInitialized!(fg)


global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



# X1 = getVal(fg, :x1)
@show global x0 = getKDEMax(getVertKDE(fg, :x0))

@test norm(x0[1:2] - [0.0;0.0]) < 1.0
@test norm(x0[3:4] - [10.0;10.0]) < 1.0


@show global x1 = getKDEMax(getVertKDE(fg, :x1))

@test norm(x1[1:2] - [20.0;20.0]) < 1.0
@test norm(x1[3:4] - [10.0;10.0]) < 1.0


end



@testset "test DynPoint2VelocityPrior..." begin

# N = 75
global fg = initfg()

# add two point locations
global v0 = addNode!(fg, :x0, DynPoint2(ut=0))
global v1 = addNode!(fg, :x1, DynPoint2(ut=1000_000))
global v2 = addNode!(fg, :x2, DynPoint2(ut=2000_000))


# Prior factor as boundary condition
global pp0 = DynPoint2VelocityPrior(MvNormal([zeros(2);10*ones(2)], 0.1*Matrix{Float64}(LinearAlgebra.I, 4,4)))
global f0 = addFactor!(fg, [:x0;], pp0)

# conditional likelihood between Dynamic Point2
global dp2dp2 = VelPoint2VelPoint2(MvNormal([10*ones(2);zeros(2)], 0.1*Matrix{Float64}(LinearAlgebra.I, 4,4)))
global f1 = addFactor!(fg, [:x0;:x1], dp2dp2)

# conditional likelihood between Dynamic Point2
global dp2dp2 = VelPoint2VelPoint2(MvNormal([10*ones(2);zeros(2)], 0.1*Matrix{Float64}(LinearAlgebra.I, 4,4)))
global f2 = addFactor!(fg, [:x1;:x2], dp2dp2)


# Graphs.plot(fg.g)
ensureAllInitialized!(fg)

@test isInitialized(fg, :x0)
@test isInitialized(fg, :x1)
@test isInitialized(fg, :x2)


global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



# X1 = getVal(fg, :x1)
@show global x0 = getKDEMax(getVertKDE(fg, :x0))
@test norm(x0[1:2] - [0.0;0.0]) < 1.0
@test norm(x0[3:4] - [10.0;10.0]) < 1.0

@show global x1 = getKDEMax(getVertKDE(fg, :x1))
@test norm(x1[1:2] - [10.0;10.0]) < 1.0
@test norm(x1[3:4] - [10.0;10.0]) < 1.0

@show global x2 = getKDEMax(getVertKDE(fg, :x2))
@test norm(x2[1:2] - [20.0;20.0]) < 1.0
@test norm(x2[3:4] - [10.0;10.0]) < 1.0


end






@testset "test VelPoint2VelPoint2" begin


global N=100
global pμ = [0.0,0,10.0,0]
global pσ = Matrix(Diagonal([0.1;0.1;0.1;0.1].^2))

global fg = initfg();

addNode!(fg, :x1, DynPoint2(ut=0))
global pp = DynPoint2VelocityPrior(MvNormal(pμ,pσ))
addFactor!(fg, [:x1;], pp, autoinit=false)

ensureAllInitialized!(fg)

addNode!(fg, :x2, DynPoint2(ut=1_000_000))
global dpμ = [10.0;0;0;0];
global dpσ = Matrix(Diagonal([1.0;1;0.5;0.01].^2))
global vp = VelPoint2VelPoint2(MvNormal(dpμ,dpσ))
addFactor!(fg, [:x1,:x2], vp, autoinit=false)

# writeGraphPdf(fg)

# note, nothign is initialized yet...
@test isInitialized(fg, :x1)
@test !isInitialized(fg, :x2)

# lets init the first
# IncrementalInference.doautoinit!(fg, :x1)

global X1 = getVal(fg, :x1)
@test 0.7*N < sum(-0.5 .< X1[1,:] .< 0.5)
@test 0.7*N < sum(-0.5 .< X1[2,:] .< 0.5)
@test 0.7*N < sum(9.0 .< X1[3,:] .< 11.0)
@test 0.7*N < sum(9.0 .< X1[3,:] .< 11.0)

# now test "projection" / probabilistic-convolution through VelPointVelPoint{D, N}
global pX2 = approxConv(fg, :x1x2f1, :x2)

@test 0.7*N < sum(9 .< pX2[1,:] .< 11)
@test 0.7*N < sum(-2 .< pX2[2,:] .< 2)
@test 0.7*N < sum(9.0 .< pX2[3,:] .< 11.0)
@test 0.7*N < sum(9.0 .< pX2[3,:] .< 11.0)

# using RoMEPlotting, Gadfly
# using KernelDensityEstimate
# using KernelDensityEstimatePlotting
# X1 = getVal(fg,:x1)
# plot(x=pX2[1,:], y=pX2[2,:], Geom.hexbin)
# pl = plotKDE(kde!(pX2[1:2,:]), dims=[1;2], levels=3, c=["blue"])
# pl = plotKDE(kde!(X1[1:2,:]), dims=[1;2], levels=3, c=["blue"])

end






#
