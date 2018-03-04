
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


N = 100
fg = initfg()

initCov = eye(6)
[initCov[i,i] = 0.01^2 for i in 4:6];
odoCov = deepcopy(initCov)


println("Adding PriorPose3 to graph...")
v1 = addNode!(fg,:x1,  0.1*randn(6,N),  N=N)
initPosePrior = PriorPose3( MvNormal( zeros(6), initCov) )
f1  = addFactor!(fg,[v1], initPosePrior)

println("Ensure vertex initialized properly")
# start with to tight an initialization
muX1 = Base.mean(getVal(fg,:x1),2)
stdX1 = Base.std(getVal(fg,:x1),2)
@test sum(map(Int,abs.(muX1) .< 0.1)) == 6
@test sum(map(Int, 0.05 .< stdX1 .< 0.15)) == 6


println("Testing PriorPose3 evaluation...")
priorpts = evalFactor2(fg, fg.g.vertices[2], 1)
means = Base.mean(priorpts,2)
@test sum(map(Int,abs.(means[1:3]) .> 0.5)) == 0
@test sum(map(Int,abs.(means[4:6]) .> 0.05)) == 0

v2, f2 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )
v3, f3 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

println("Testing Pose3Pose3 evaluation...")
X1pts = getVal(fg, :x1)
X2pts = evalFactor2(fg, fg.g.vertices[6], 3, N=N)
X3pts = evalFactor2(fg, fg.g.vertices[6], 5)
X2ptsMean = Base.mean(X2pts,2)
X3ptsMean = Base.mean(X3pts,2)

@test  sum(map(Int, abs.(X2ptsMean) - [25.0;0;0;0;0;0] .< 5.0 )) == 6
@test  sum(map(Int, abs.(X3ptsMean - [50.0;0;0;0;0;0]) .< 5.0 )) == 6

tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg,tree,N=N)
# inferOverTree!(fg,tree)

println("Adding Pose3Pose3NH to graph...")

odo3 = SE3([-35;0;0], Quaternion(0))
odoc3 = Pose3Pose3NH( MvNormal(veeEuler(odo3), odoCov), [0.5;0.5]) # define 50/50% hypothesis
addFactor!(fg,[v3;v1],odoc3)


X1pts = evalFactor2(fg, fg.g.vertices[7], 1, N=N)
X2pts = evalFactor2(fg, fg.g.vertices[7], 5, N=N)

p1 = kde!(X1pts)
p2 = kde!(X2pts)


# plotKDE([marginal(p1,[1]),marginal(p1t,[1])], c=["black";"red"])


using JLD, HDF5

println("loading validation data for testing.")
@load joinpath(dirname(@__FILE__),"testvalidation.jld") X1ptst X2ptst
# @save joinpath(dirname(@__FILE__),"testvalidation.jld") X1ptst X2ptst

p1t = kde!(X1ptst)
p2t = kde!(X2ptst)

# plotKDE([p2t;p2],c=["red";"blue"],dims=[1;2],levels=3)
# kld(marginal(p1,[2]), marginal(p1t,[2]), method=:unscented)

t1 = minimum([abs(kld(p1, p1t)[1]) ; abs(kld(p1t, p1)[1])])
t2 = minimum([abs(kld(p2, p2t)[1]) ; abs(kld(p2t, p2)[1])])

@test t1 < 30.0
@test t2 < 30.0


# plotKDE(margisal(p1,[1]))

# warn("incomplete Pose3Pose3NH tests.")
