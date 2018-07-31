using RoME, IncrementalInference, Distributions
using Base: Test


@testset "test DynPose2 and velocity..." begin

N = 75
fg = initfg()

# add first pose locations
addNode!(fg, :x0, DynPose2(ut=0))

# Prior factor as boundary condition
pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), diagm([0.01; 0.01; 0.001].^2)),
                            MvNormal(zeros(2), diagm([0.1; 0.1].^2)))
addFactor!(fg, [:x0;], pp0)

# initialize the first pose
doautoinit!(fg, :x0)

addNode!(fg, :x1, DynPose2(ut=1000_000))

# conditional likelihood between Dynamic Point2
dp2dp2 = VelPose2VelPose2(MvNormal([10.0;0;0], diagm([0.01;0.01;0.001].^2)),
                          MvNormal(zeros(2), diagm([0.1; 0.1].^2)))
addFactor!(fg, [:x0;:x1], dp2dp2)

pts = approxConv(fg, :x0x1f1, :x1)

# Graphs.plot(fg.g)
ensureAllInitialized!(fg)




tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)

end








#
