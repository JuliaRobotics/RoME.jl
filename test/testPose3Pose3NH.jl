
using RoME
# , Distributions
# using KernelDensityEstimate
using Test
using DelimitedFiles
using Statistics

@testset "Test syntax for null hypothesis matrix substitution" begin
  function ff(X::Array{Float64,2})
    return X.+1.0
  end

  global c = Categorical([0.5;0.5])
  global A = zeros(3,5)

  global skipdos = rand(c, 5)
  global dos = skipdos.==2

  global B = zeros(3,5)
  B[:,dos] = ff(A[:,dos])

  @test sum(sum(B[:,dos],dims=2) .== sum(dos),dims=1)[1,1] == 3
end


println("Test if null hypothesis occurs as expected...")


global N = 100
global fg = initfg()

global initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.01^2 for i in 4:6];
global odoCov = deepcopy(initCov)


println("Adding PriorPose3 to graph...")
global v1 = addNode!(fg, :x1, Pose3, N=N)
global initPosePrior = PriorPose3( MvNormal( zeros(6), initCov) )
global f1  = addFactor!(fg,[v1], initPosePrior, autoinit=true)

# ls(fg, :x1)

ensureAllInitialized!(fg)

println("Ensure vertex initialized properly")
# start with to tight an initialization
global muX1 = Statistics.mean(getVal(fg,:x1),dims=2)
global stdX1 = Statistics.std(getVal(fg,:x1),dims=2)
@test sum(map(Int,abs.(muX1[1:3]) .< 0.4)) == 3
@test sum(map(Int,abs.(muX1[4:6]) .< 0.04)) == 3
@test sum(map(Int,abs.(1.0 .- stdX1[1:3]) .< 0.3)) == 3
@test sum(map(Int,abs.(0.01 .- stdX1[4:6]) .< 0.1)) == 3


println("Testing PriorPose3 evaluation...")
global priorpts = approxConv(fg, :x1f1, :x1)
# priorpts = evalFactor2(fg, fg.g.vertices[2], 1)
global means = Statistics.mean(priorpts,dims=2)
@test sum(map(Int,abs.(means[1:3]) .> 0.5)) == 0
@test sum(map(Int,abs.(means[4:6]) .> 0.05)) == 0

global v2, f2 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )
global v3, f3 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

println("Testing Pose3Pose3 evaluation...")
global X1pts = getVal(fg, :x1)
global X2pts = approxConv(fg, :x2x3f1, :x2, N=N)
# X2pts = evalFactor2(fg, fg.g.vertices[6], 3, N=N)
global X3pts = approxConv(fg, :x2x3f1, :x3)
global X2ptsMean = Statistics.mean(X2pts,dims=2)
global X3ptsMean = Statistics.mean(X3pts,dims=2)

@test  sum(map(Int, abs.(X2ptsMean) - [25.0;0;0;0;0;0] .< 5.0 ))  == 6
@test  sum(map(Int, abs.(X3ptsMean -  [50.0;0;0;0;0;0]) .< 5.0 )) == 6


global tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg,tree,N=N)
# inferOverTree!(fg,tree)

println("Adding Pose3Pose3NH to graph...")

global odo3 = SE3([-35;0;0], Quaternion(0))
global odoc3 = Pose3Pose3NH( MvNormal(veeEuler(odo3), odoCov), [0.5;0.5]) # define 50/50% hypothesis

@test ls(fg)[1] == [:x1;:x2;:x3]
addFactor!(fg,[:x3;:x1],odoc3)

global X1pts = approxConv(fg, :x3x1f1, :x1, N=N)
global X2pts = approxConv(fg, :x3x1f1, :x3, N=N)

global p1 = kde!(X1pts)
global p2 = kde!(X2pts)


@testset "loading validation data for testing." begin
    # @save joinpath(dirname(@__FILE__),"testvalidation.jld") X1ptst X2ptst
    global tstdtdir = dirname(@__FILE__)
    global X1ptst = readdlm(joinpath(tstdtdir, "X1ptst.csv"),',')
    global X2ptst = readdlm(joinpath(tstdtdir, "X2ptst.csv"),',')

    global p1t = kde!(X1ptst)
    global p2t = kde!(X2ptst)

    # plotKDE([p2t;p2],c=["red";"blue"],dims=[1;2],levels=3)
    # kld(marginal(p1,[2]), marginal(p1t,[2]), method=:unscented)

    global t1 = minimum([abs(kld(p1, p1t)[1]) ; abs(kld(p1t, p1)[1])])
    global t2 = minimum([abs(kld(p2, p2t)[1]) ; abs(kld(p2t, p2)[1])])

    @test t1 < 80.0
    @test t2 < 80.0
end



# using KernelDensityEstimatePlotting
# plotKDE([p1; p1t], dims=[1], c=["black";"red"])
