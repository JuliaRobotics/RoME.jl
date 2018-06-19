using RoME, IncrementalInference, TransformUtils
using Distributions
using KernelDensityEstimate
using Base.Test


tf = SE3([0.0;0.0;0.0], AngleAxis(pi/4,[0;0;1.0]))# Euler(pi/4,0.0,0.0) )

N = 1
initCov = 0.0001*eye(6)
[initCov[i,i] = 0.000001 for i in 4:6];
odoCov = deepcopy(initCov)
odo = Pose3Pose3( MvNormal(veeEuler(tf), odoCov) )

X = [0.01*randn(5,N);0*pi/4+0.01*randn(1,N)]

Y = X ⊕ odo

@test norm(Y[1:3]-zeros(3)) < 1.0
@test norm(Y[4:6]-[zeros(2);pi/4]) < 0.15



tf = SE3([0.0;0.0;0.0], AngleAxis(pi/4,[0;0;1.0]))# Euler(pi/4,0.0,0.0) )

N = 1
initCov = 0.01*eye(6)
[initCov[i,i] = 0.001 for i in 4:6];
odoCov = deepcopy(initCov)
odo = Pose3Pose3(  MvNormal(veeEuler(tf), odoCov) )

X = [0.01*randn(5,N);0*pi/4+0.01*randn(1,N)]

Y = X ⊕ odo

@test norm(Y[1:3]-zeros(3)) < 1.0
@test norm(Y[4:6]-[zeros(2);pi/4]) < 0.2




N = 75
fg = initfg()

initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)


@testset "Adding PriorPose3 to graph..." begin
  v1 = addNode!(fg, :x1, Pose3,  N=N) # 0.1*randn(6,N)
  initPosePrior = PriorPose3( MvNormal(zeros(6), initCov) )
  f1  = addFactor!(fg,[:x1;], initPosePrior)
  @test !isInitialized(fg, :x1)
end


@testset "Ensure vertex initialized properly" begin
  # start with initialization
  ensureAllInitialized!(fg)
  @test isInitialized(fg, :x1)
  @show muX1 = Base.mean(getVal(fg,:x1),2)
  @show stdX1 = Base.std(getVal(fg,:x1),2)
  @test sum(map(Int,abs.(muX1[1:3]) .< 0.5)) == 3
  @test sum(map(Int,abs.(muX1[4:6]) .< 0.05)) == 3
  @test sum(map(Int, 0.5 .< stdX1[1:3] .< 1.5)) == 3
  @test sum(map(Int, 0.05 .< stdX1[4:6] .< 0.15)) == 3
end


@testset "Testing PriorPose3 evaluation..." begin
  priorpts = evalFactor2(fg, fg.g.vertices[2], 1)
  means = Base.mean(priorpts,2)
  @test sum(map(Int,abs.(means[1:3]) .> 0.5)) == 0
  @test sum(map(Int,abs.(means[4:6]) .> 0.05)) == 0
end


@testset "Adding Pose3Pose3 to graph..." begin
  odo = SE3([10;0;0], Quaternion(0))
  pts0X2 = projectParticles(getVal(fg,:x1), MvNormal(veeEuler(odo), odoCov) )
  odoconstr = Pose3Pose3( MvNormal(veeEuler(odo), odoCov) )
  v2 = addNode!(fg,:x2, Pose3, N=N) # pts0X2
  addFactor!(fg,[:x1;:x2],odoconstr)
  @test !isInitialized(fg, :x2)
end


@testset "Testing Pose3Pose3 evaluation..." begin
  ensureAllInitialized!(fg)
  @test isInitialized(fg, :x2)
  X1pts = evalFactor2(fg, fg.g.vertices[4], 1)
  X2pts = evalFactor2(fg, fg.g.vertices[4], 3)
  X2ptsMean = Base.mean(X2pts,2)
  X1ptsMean = Base.mean(X1pts,2)
  @show X1ptsMean
  @test  sum(map(Int, abs.(X1ptsMean) .< 1.0 )) == 6
  @test  sum(map(Int, abs.(X2ptsMean - [10.0;0;0;0;0;0]) .< 1.0 )) == 6
end


@testset "Construct Bayes tree and perform inference..." begin
  tree = prepBatchTree!(fg);
  inferOverTree!(fg, tree, N=N)
  @test true
end

@testset "Ensure basic parameters on x1,x2 after inference..." begin
  # check mean and covariances after one up and down pass over the tree
  muX1 = Base.mean(getVal(fg,:x1),2)
  stdX1 = Base.std(getVal(fg,:x1),2)
  @test sum(map(Int,abs.(muX1[1:3]) .< 1.0)) == 3
  @test sum(map(Int,abs.(muX1[4:6]) .< 0.1)) == 3
  @test sum(map(Int, 0.4 .< stdX1[1:3] .< 1.6)) == 3 # had a 2==3 failure here
  @test sum(map(Int, 0.025 .< stdX1[4:6] .< 0.25)) == 3
  muX2 = Base.mean(getVal(fg,:x2),2)
  stdX2 = Base.std(getVal(fg,:x2),2)
  @show muX2[1:3]-[10.0;0;0]
  @test sum(map(Int, abs.(muX2[1:3]-[10.0;0;0]) .< 1.5)) == 3
  @test sum(map(Int, abs.(muX2[4:6]) .< 0.1)) == 3
  @show println("previous test failure 0.75 .< $(round.(stdX2[1:3],2)) .< 2.25")
  @test sum(map(Int, 0.75 .< stdX2[1:3] .< 2.25)) == 3
  @test sum(map(Int, 0.05 .< stdX2[4:6] .< 0.25)) == 3
end

# println("Plot marginals to see what is happening")
# plotKDE(marginal(getVertKDE(fg,:x1),[1]))
# plotKDE(marginal(getVertKDE(fg,:x2),[1]))




#
