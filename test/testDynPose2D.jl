using RoME, IncrementalInference, Distributions
using Base: Test

const TU = TransformUtils


@testset "test DynPose2 and velocity..." begin

N = 75
fg = initfg()

# add first pose locations
addNode!(fg, :x0, DynPose2(ut=0))

# Prior factor as boundary condition
pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), diagm([0.01; 0.01; 0.001].^2)),
                            MvNormal([10.0;0], diagm([0.1; 0.1].^2)))
addFactor!(fg, [:x0;], pp0)

# initialize the first pose
IncrementalInference.doautoinit!(fg, [getVert(fg,:x0);])

addNode!(fg, :x1, DynPose2(ut=1000_000))

# conditional likelihood between Dynamic Point2
dp2dp2 = VelPose2VelPose2(MvNormal([10.0;0;0], diagm([0.01;0.01;0.001].^2)),
                          MvNormal([0.0;0], diagm([0.1; 0.1].^2)))
addFactor!(fg, [:x0;:x1], dp2dp2)

# getVal(fg,:x0)
pts = approxConv(fg, :x0x1f1, :x1)

# Graphs.plot(fg.g)
# ensureAllInitialized!(fg)

tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)

X1 = getVal(fg, :x1)

@test 0.9*N <= sum(abs.(X1[1,:] - 10.0) .< 0.5)
@test 0.9*N <= sum(abs.(X1[2,:] - 0.0) .< 0.5)
@show TU.wrapRad.(X1[3,:])
@test 0.8*N <= sum(abs.(TU.wrapRad.(X1[3,:]) - 0.0) .< 0.1)
warn("wrapRad issue, accepting 80% as good enough until issue JuliaRobotics/RoME.jl#90 is fixed.")
@test 0.9*N <= sum(abs.(X1[4,:] - 10.0) .< 0.5)
@test 0.9*N <= sum(abs.(X1[5,:] - 0.0) .< 0.5)


end


@testset "test distribution compare functions..." begin

mu = randn(6)
mv1 = MvNormal(deepcopy(mu), eye(6))
mv2 = MvNormal(deepcopy(mu), eye(6))
mv3 = MvNormal(randn(6), eye(6))
@test RoME.compare(mv1, mv2)
@test !RoME.compare(mv1, mv3)
@test !RoME.compare(mv2, mv3)

end


@testset "test DynPose2 packing converters..." begin

pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), diagm([0.01; 0.01; 0.001].^2)),
                            MvNormal([10.0;0], diagm([0.1; 0.1].^2)))

pp = convert(PackedDynPose2VelocityPrior, pp0)
ppu = convert(DynPose2VelocityPrior, pp)

@test RoME.compare(pp0, ppu)

dp2dp2 = VelPose2VelPose2(MvNormal([10.0;0;0], diagm([0.01;0.01;0.001].^2)),
                          MvNormal([0.0;0], diagm([0.1; 0.1].^2)))

pp = convert(PackedVelPose2VelPose2, dp2dp2)
ppu = convert(VelPose2VelPose2, pp)

@test RoME.compare(dp2dp2, ppu)

end






#
