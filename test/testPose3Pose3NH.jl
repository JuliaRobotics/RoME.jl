
using RoME, IncrementalInference, TransformUtils, Distributions
using KernelDensityEstimate
using Base.Test


begin
  function ff(X::Array{Float64,2})
    return X.+1.0
  end

  c = Categorical([0.5;0.5])
  A = zeros(3,5)

  skipdos = rand(c, 5)
  dos = skipdos.==2

  B = zeros(3,5)
  B[:,dos] = ff(A[:,dos])

  @test sum(sum(B[:,dos],2) .== sum(dos),1)[1,1] == 3
  println("Syntax for null hypothesis matrix substitutions work.")
end


println("Test if null hypothesis occurs as expected...")


N = 200
fg = initfg()

initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)


println("Adding PriorPose3 to graph...")
v1 = addNode!(fg,"x1",  0.1*randn(6,N),  N=N)
initPosePrior = PriorPose3(SE3(0), initCov)
f1  = addFactor!(fg,[v1], initPosePrior)

println("Ensure vertex initialized properly")
# start with to tight an initialization
muX1 = Base.mean(getVal(fg,"x1"),2)
stdX1 = Base.std(getVal(fg,"x1"),2)
@test sum(map(Int,abs(muX1) .< 0.1)) == 6
@test sum(map(Int, 0.05 .< stdX1 .< 0.15)) == 6


println("Testing PriorPose3 evaluation...")
priorpts = evalFactor2(fg, fg.g.vertices[2], 1)
means = Base.mean(priorpts,2)
@test sum(map(Int,abs(means[1:3]) .> 0.5)) == 0
@test sum(map(Int,abs(means[4:6]) .> 0.05)) == 0



println("Adding Pose3Pose3NH to graph...")
odo = SE3([10;0;0], Quaternion(0))
pts0X2 = projectParticles(getVal(fg,"x1"), odo, odoCov)
odoconstr = Pose3Pose3NH(odo, odoCov, [0.5;0.5]) # define 50/50% hypothesis
v2 = addNode!(fg,"x2",  pts0X2, N=N)
addFactor!(fg,[v1;v2],odoconstr)



println("Testing Pose3Pose3NH evaluation...")
X1pts = evalFactor2(fg, fg.g.vertices[4], 1)
X2pts = evalFactor2(fg, fg.g.vertices[4], 3)
# X2ptsMean = Base.mean(X2pts,2)
# X1ptsMean = Base.mean(X1pts,2)
# @test  sum(map(Int, abs(X1ptsMean) .< 0.5 )) == 6
# @test  sum(map(Int, abs(X2ptsMean - [10.0;0;0;0;0;0]) .< 0.5 )) == 6




# warn("incomplete Pose3Pose3NH tests.")
