using RoME
using Test
# using RoMEPlotting, Distributions

import IncrementalInference: getSample

mutable struct NorthSouthPartial{T} <: FunctorSingleton
  Z::T
  partial::Tuple
  NorthSouthPartial{D}() where D = new{D}()
  NorthSouthPartial{D}(Z::D) where {D <: IIF.SamplableBelief} = new{D}(Z, (2,))
end
NorthSouthPartial(Z::D) where {D <: IIF.SamplableBelief} = NorthSouthPartial{D}(Z)

getSample(ns::NorthSouthPartial, N=1) = (reshape(rand(ns.Z, N),1,N),)



@testset "test standard multimodal range bearing factor setup..." begin

# Start with an empty graph
global N = 100
global fg = initfg(sessionname="MULTIMODAL_2D_TUTORIAL")

# Add landmarks with Bearing range measurements
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
addNode!(fg, :l2, Point2, labels=["LANDMARK"])

addFactor!(fg, [:l1], Prior(MvNormal([10.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) )
addFactor!(fg, [:l2], Prior(MvNormal([30.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) )

setVal!(fg, :l1, predictbelief(fg, :l1, [:l1f1;]))
setVal!(fg, :l2, predictbelief(fg, :l2, [:l2f1;]))

addNode!(fg, :x0, Pose2)
# addFactor!(fg, [:x0], Prior(MvNormal([0.0;0.0;0], Matrix(Diagonal([1.0;1.0;0.01].^2)))) )
addFactor!(fg, [:x0;], NorthSouthPartial(Normal(0,1.0)))


global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1; :l2], p2br, multihypo=[1.0; 0.5; 0.5])

predictbelief(fg, :x0, ls(fg, :x0))




@test true
# writeGraphPdf(fg)
end





@testset "test multimodal bearing range factors calculate pose position properly..." begin

ensureAllInitialized!(fg)

global tree = wipeBuildNewTree!(fg)

global infv = [inferOverTreeR!(fg, tree, N=N) for i in 1:4]
@test length(infv) == 4


# X0 = getVertKDE(fg, :x0)
# X0pts = getPoints(X0)
global X0pts = getVal(fg, :x0)


@test size(X0pts,1) == 3
@test size(X0pts,2) == N

@test sum([0 < sum(-20 .< X0pts[1,:] .< 0);
           0 < sum(  0 .< X0pts[1,:] .< 20);
           0 < sum( 20 .< X0pts[1,:] .< 40);
           0 < sum( 40 .< X0pts[1,:] .< 60)]) > 2

end


@testset "test multimodal landmark locations are computed correclty..." begin


# Start with an empty graph
global fg = initfg(sessionname="MULTIMODAL_2D_TUTORIAL")

addNode!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal([0.0;0.0;0], Matrix(Diagonal([1.0;1.0;0.01].^2)))) ) # TODO IIF.Prior with IIF 0.3.9

# Add landmarks with Bearing range measurements
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
addFactor!(fg, [:l1], PriorPose2(MvNormal([40.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) ) # TODO IIF.Prior with IIF 0.3.9

addNode!(fg, :l2, Point2, labels=["LANDMARK"])
addFactor!(fg, [:l2;], NorthSouthPartial(Normal(0,1.0)))
# addFactor!(fg, [:l2], PriorPose2(MvNormal([30.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) ) # TODO IIF.Prior with IIF 0.3.9

global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1; :l2], p2br, multihypo=[1.0; 0.5; 0.5])

@test true

# setVal!(getVert(fg, :l2), zeros(2,N))
# writeGraphPdf(fg)
global tree = wipeBuildNewTree!(fg)

global infv = [inferOverTreeR!(fg, tree, N=N) for i in 1:4]
@test length(infv) == 4


# L1 = getVertKDE(fg, :l1)
global L2 = getVertKDE(fg, :l2)
# X0 = getVertKDE(fg, :x0)

global L2pts = getPoints(L2)

@test size(L2pts,1) == 2
@test size(L2pts,2) == N

# plot(L2, dims=[1;2])

global numM1 = sum(10 .< L2pts[1,:] .< 30)
global numM2 = sum(30 .< L2pts[1,:] .< 50)

@test 30 < numM1 < 70
@test 30 < numM2 < 70
@test numM1 + numM2 == N

end
