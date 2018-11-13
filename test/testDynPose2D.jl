using RoME
# , IncrementalInference, Distributions
using Test

# const TU = TransformUtils


@testset "test DynPose2 and velocity..." begin

global N = 100
global fg = initfg()

# add first pose locations
addNode!(fg, :x0, DynPose2(ut=0))

# Prior factor as boundary condition
global pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), Matrix(Diagonal([0.01; 0.01; 0.001].^2))),
                            MvNormal([10.0;0], Matrix(Diagonal([0.1; 0.1].^2))))
addFactor!(fg, [:x0;], pp0)

# initialize the first pose
IncrementalInference.doautoinit!(fg, [getVert(fg,:x0);])

addNode!(fg, :x1, DynPose2(ut=1000_000))

# conditional likelihood between Dynamic Point2
global dp2dp2 = VelPose2VelPose2(MvNormal([10.0;0;0], Matrix(Diagonal([0.01;0.01;0.001].^2))),
                          MvNormal([0.0;0], Matrix(Diagonal([0.1; 0.1].^2))))
addFactor!(fg, [:x0;:x1], dp2dp2)

# getVal(fg,:x0)
global pts = approxConv(fg, :x0x1f1, :x1)

# Graphs.plot(fg.g)
# ensureAllInitialized!(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)

global X1 = getVal(fg, :x1)

@test 0.9*N <= sum(abs.(X1[1,:] .- 10.0) .< 0.5)
@test 0.9*N <= sum(abs.(X1[2,:] .- 0.0) .< 0.5)
@show TU.wrapRad.(X1[3,:])
@test 0.8*N <= sum(abs.(TU.wrapRad.(X1[3,:]) .- 0.0) .< 0.1)
@warn "wrapRad issue, accepting 80% as good enough until issue JuliaRobotics/RoME.jl#90 is fixed."
@test 0.9*N <= sum(abs.(X1[4,:] .- 10.0) .< 0.5)
@test 0.9*N <= sum(abs.(X1[5,:] .- 0.0) .< 0.5)


end


@testset "test distribution compare functions..." begin

global mu = randn(6)
global mv1 = MvNormal(deepcopy(mu), Matrix{Float64}(LinearAlgebra.I, 6,6))
global mv2 = MvNormal(deepcopy(mu), Matrix{Float64}(LinearAlgebra.I, 6,6))
global mv3 = MvNormal(randn(6), Matrix{Float64}(LinearAlgebra.I, 6,6))
@test RoME.compare(mv1, mv2)
@test !RoME.compare(mv1, mv3)
@test !RoME.compare(mv2, mv3)

end


@testset "test DynPose2 packing converters..." begin

global pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), Matrix(Diagonal([0.01; 0.01; 0.001].^2))),
                            MvNormal([10.0;0], Matrix(Diagonal([0.1; 0.1].^2))))

global pp = convert(PackedDynPose2VelocityPrior, pp0)
global ppu = convert(DynPose2VelocityPrior, pp)

@test RoME.compare(pp0, ppu)

global dp2dp2 = VelPose2VelPose2(MvNormal([10.0;0;0], Matrix(Diagonal([0.01;0.01;0.001].^2))),
                          MvNormal([0.0;0], Matrix(Diagonal([0.1; 0.1].^2))))

global pp = convert(PackedVelPose2VelPose2, dp2dp2)
global ppu = convert(VelPose2VelPose2, pp)

@test RoME.compare(dp2dp2, ppu)

end



@testset "test many DynPose2 chain stationary and 'pulled'..." begin

global N = 100
global fg = initfg()

# add first pose locations
addNode!(fg, :x0, DynPose2(ut=0))

# Prior factor as boundary condition
global pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), Matrix(Diagonal([0.01; 0.01; 0.001].^2))),
                            MvNormal([0.0;0], Matrix(Diagonal([0.1; 0.1].^2))))
addFactor!(fg, [:x0;], pp0)

global sym = :x0
global k = 0
for sy in Symbol[Symbol("x$i") for i in 1:10]

global k+=1
addNode!(fg, sy, DynPose2(ut=1000_000*k))

# conditional likelihood between Dynamic Point2
global dp2dp2 = VelPose2VelPose2(MvNormal([0.0;0;0], Matrix(Diagonal([1.0;0.1;0.001].^2))),
                          MvNormal([0.0;0], Matrix(Diagonal([0.1; 0.1].^2))))
addFactor!(fg, [sym;sy], dp2dp2)
global sym =sy

end


global x5 = KDE.getKDEMean(getVertKDE(fg, :x5))

@test abs(x5[1]) < 1.25
@test abs(x5[2]) < 1.25
@test abs(TU.wrapRad(x5[3])) < 0.4
@test abs(x5[4]) < 0.5
@test abs(x5[5]) < 0.5


ensureAllInitialized!(fg)

global x10 = KDE.getKDEMean(getVertKDE(fg, :x10))

@test abs(x10[1]) < 1.25
@test abs(x10[2]) < 1.25
@test abs(TU.wrapRad(x10[3])) < 0.4
@test abs(x10[4]) < 0.5
@test abs(x10[5]) < 0.5

# using RoMEPlotting
# drawPoses(fg)
# plotPose(fg, [:x10])

# batchSolveR!(fg, N=N)
tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree, N=N)


global x5 = KDE.getKDEMean(getVertKDE(fg, :x5))

@test abs(x5[1]) < 1.5
@test abs(x5[2]) < 1.5
@test abs(TU.wrapRad(x5[3])) < 0.4
@test abs(x5[4]) < 0.5
@test abs(x5[5]) < 0.5

global x10 = KDE.getKDEMean(getVertKDE(fg, :x10))

@test abs(x10[1]) < 2.75
@test abs(x10[2]) < 2.75
@test abs(TU.wrapRad(x10[3])) < 0.5
@test abs(x10[4]) < 0.5
@test abs(x10[5]) < 0.5



# pull the tail end out with position
global pp10 = DynPose2VelocityPrior(MvNormal([10.0;0;0], Matrix(Diagonal([0.01; 0.01; 0.001].^2))),
                            MvNormal([0.0;0], Matrix(Diagonal([0.1; 0.1].^2))))
addFactor!(fg, [:x10;], pp10)



batchSolve!(fg, N=N)
# run(`evince /tmp/bt.pdf`)



global x10 = KDE.getKDEMean(getVertKDE(fg, :x10))

@test 5.0 < x10[1]
@test abs(x10[2]) < 1.0
@test abs(TU.wrapRad(x10[3])) < 0.6
@test -0.1 < x10[4] < 1.0
@test abs(x10[5]) < 0.5


for sym in [Symbol("x$i") for i in 2:9]

global XX = KDE.getKDEMean(getVertKDE(fg, sym))

@show sym, round.(XX,digits=5)
@test -1.5 < XX[1] < 10.0
@test abs(XX[2]) < 1.0
@test abs(TU.wrapRad(XX[3])) < 1.0
@test -0.3 < XX[4] < 2.0
@test abs(XX[5]) < 0.5

end

end


# using RoMEPlotting
# plotLocalProduct(fg, :x10, dims=[1;2])
# drawPoses(fg)
# plotPose(fg, [:x9],levels=1);
#
# plotPose(fg, [:x1;:x2;:x3;:x4;:x5;:x6;:x7;:x8;:x9;:x10],levels=1);

# savejld(fg) # tempfg.jld



@testset "test many DynPose2 sideways velocity..." begin

global N = 100
global fg = initfg()

# add first pose locations
addNode!(fg, :x0, DynPose2(ut=0))

# Prior factor as boundary condition
global pp0 = DynPose2VelocityPrior(MvNormal([0.0;0.0;pi/2], Matrix(Diagonal([0.01; 0.01; 0.001].^2))),
                            MvNormal([0.0;0], Matrix(Diagonal([0.5; 0.5].^2))))
addFactor!(fg, [:x0;], pp0)


addNode!(fg, :x1, DynPose2(ut=1000_000))

global pp0 = DynPose2VelocityPrior(MvNormal([1.0;0.0;pi/2], Matrix(Diagonal([0.01; 0.01; 0.001].^2))),
                            MvNormal([0.0;0], Matrix(Diagonal([0.5; 0.5].^2))))
addFactor!(fg, [:x1;], pp0)

# conditional likelihood between Dynamic Point2
global dp2dp2 = VelPose2VelPose2(MvNormal([0.0;-1.0;0], Matrix(Diagonal([0.01;0.01;0.001].^2))),
                          MvNormal([0.0;0], Matrix(Diagonal([0.1; 0.1].^2))))
addFactor!(fg, [:x0;:x1], dp2dp2)


batchSolve!(fg,N=N)


# test for velocity in the body frame
global x0 = KDE.getKDEMean(getVertKDE(fg, :x0))

@test -0.4 < x0[1] < 2.0
@test abs(x0[2]) < 0.5
@test abs(x0[3] - pi/2) < 0.1
@test abs(x0[4]) < 0.4
@test -1.5 < x0[5] < -0.5


global x1 = KDE.getKDEMean(getVertKDE(fg, :x1))

@test -0.1 < x1[1] < 2.0
@test abs(x1[2]) < 0.5
@test abs(x1[3] - pi/2) < 0.1
@test abs(x1[4]) < 0.4
@test -1.5 < x1[5] < -0.5


end

# using RoMEPlotting
#
# drawPoses(fg)
#
# plotPose(fg, [:x0;:x1]);






#
