using RoME
using TransformUtils
using Statistics
# , IncrementalInference, TransformUtils
# using Distributions
# using KernelDensityEstimate
using Test

# const TU = TransformUtils

global tf = SE3([0.0;0.0;0.0], TU.AngleAxis(pi/4,[0;0;1.0]))# Euler(pi/4,0.0,0.0) )

global N = 1
global initCov = 0.0001*Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.000001 for i in 4:6];
global odoCov = deepcopy(initCov)
global odo = Pose3Pose3( MvNormal(veeEuler(tf), odoCov) )

global X = [0.01*randn(5,N); (0*pi/4 .+ 0.01*randn(1,N))]

global Y = X ⊕ odo

@test norm(Y[1:3]-zeros(3)) < 1.0
@test norm(Y[4:6]-[zeros(2);pi/4]) < 0.15



global tf = SE3([0.0;0.0;0.0], TU.AngleAxis(pi/4,[0;0;1.0]))# Euler(pi/4,0.0,0.0) )

global N = 1
global initCov = 0.01*Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.001 for i in 4:6];
global odoCov = deepcopy(initCov)
global odo = Pose3Pose3(  MvNormal(veeEuler(tf), odoCov) )

global X = [0.01*randn(5,N); (0*pi/4.0 .+ 0.01*randn(1,N))]

global Y = X ⊕ odo

@test norm(Y[1:3]-zeros(3)) < 1.0
@test norm(Y[4:6]-[zeros(2);pi/4]) < 0.2




global N = 75
global fg = initfg()

global initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.01 for i in 4:6];
global odoCov = deepcopy(initCov)


@testset "Adding PriorPose3 to graph..." begin
  global v1 = addNode!(fg, :x1, Pose3,  N=N) # 0.1*randn(6,N)
  global initPosePrior = PriorPose3( MvNormal(zeros(6), initCov) )
  global f1  = addFactor!(fg,[:x1;], initPosePrior)
  @test !isInitialized(fg, :x1)
end


@testset "Ensure vertex initialized properly" begin
  # start with initialization
  ensureAllInitialized!(fg)
  @test isInitialized(fg, :x1)
  @show muX1 = Statistics.mean(getVal(fg,:x1),dims=2)
  @show stdX1 = Statistics.std(getVal(fg,:x1),dims=2)
  @test sum(map(Int,abs.(muX1[1:3]) .< 0.5)) == 3
  @test sum(map(Int,abs.(muX1[4:6]) .< 0.05)) == 3
  @test sum(map(Int, 0.5 .< stdX1[1:3] .< 1.5)) == 3
  @test sum(map(Int, 0.05 .< stdX1[4:6] .< 0.15)) == 3
end


@testset "Testing PriorPose3 evaluation..." begin
  global priorpts = evalFactor2(fg, fg.g.vertices[2], 1)
  global means = Statistics.mean(priorpts,dims=2)
  @test sum(map(Int,abs.(means[1:3]) .> 0.5)) == 0
  @test sum(map(Int,abs.(means[4:6]) .> 0.05)) == 0
end


@testset "Adding Pose3Pose3 to graph..." begin
  global odo = SE3([10;0;0], Quaternion(0))
  global pts0X2 = projectParticles(getVal(fg,:x1), MvNormal(veeEuler(odo), odoCov) )
  global odoconstr = Pose3Pose3( MvNormal(veeEuler(odo), odoCov) )
  global v2 = addNode!(fg,:x2, Pose3, N=N) # pts0X2
  addFactor!(fg,[:x1;:x2],odoconstr)
  @test !isInitialized(fg, :x2)
end


# Noticed a DomainError on convolutions here after mutlithreading upgrade.  Previously used fill(PP3REUSE, Threads.nthreads())
@testset "Testing Pose3Pose3 evaluation..." begin

ensureAllInitialized!(fg)
@test isInitialized(fg, :x2)
global X1pts = approxConv(fg, :x1x2f1, :x1)
# X1pts = evalFactor2(fg, fg.g.vertices[4], 1)
global X2pts = approxConv(fg, :x1x2f1, :x2)
# X2pts = evalFactor2(fg, fg.g.vertices[4], 3)
global X2ptsMean = Statistics.mean(X2pts,dims=2)
global X1ptsMean = Statistics.mean(X1pts,dims=2)
@show X1ptsMean
@test  sum(map(Int, abs.(X1ptsMean) .< 1.25 )) == 6
@test  sum(map(Int, abs.(X2ptsMean .- [10.0;0;0;0;0;0]) .< 1.25 )) == 6

end

@testset "Construct Bayes tree and perform inference..." begin
  global tree = wipeBuildNewTree!(fg);
  inferOverTree!(fg, tree, N=N)
  @test true
end

@testset "Ensure basic parameters on x1,x2 after inference..." begin
  # check mean and covariances after one up and down pass over the tree
  global muX1 = Statistics.mean(getVal(fg,:x1),dims=2)
  global stdX1 = Statistics.std(getVal(fg,:x1),dims=2)
  @test sum(map(Int,abs.(muX1[1:3]) .< 1.0)) == 3
  @test sum(map(Int,abs.(muX1[4:6]) .< 0.1)) == 3
  @test sum(map(Int, 0.4 .< stdX1[1:3] .< 1.6)) == 3 # had a 2==3 failure here
  @test sum(map(Int, 0.025 .< stdX1[4:6] .< 0.25)) == 3
  global muX2 = Statistics.mean(getVal(fg,:x2),dims=2)
  global stdX2 = Statistics.std(getVal(fg,:x2),dims=2)
  @show muX2[1:3]-[10.0;0;0]
  @test sum(map(Int, abs.(muX2[1:3]-[10.0;0;0]) .< 1.5)) == 3
  @test sum(map(Int, abs.(muX2[4:6]) .< 0.1)) == 3
  @show println("previous test failure 0.75 .< $(round.(stdX2[1:3],digits=2)) .< 2.25")
  @test sum(map(Int, 0.75 .< stdX2[1:3] .< 2.25)) == 3
  @test sum(map(Int, 0.05 .< stdX2[4:6] .< 0.25)) == 3
end

# println("Plot marginals to see what is happening")
# plotKDE(marginal(getVertKDE(fg,:x1),[1]))
# plotKDE(marginal(getVertKDE(fg,:x2),[1]))




#
