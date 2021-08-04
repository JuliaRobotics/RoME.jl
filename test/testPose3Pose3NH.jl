
using RoME
using Test
using DelimitedFiles
using Statistics


# TODO this initial test is likely obsolete
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


@testset "Test if null hypothesis occurs as expected..." begin

global N = 100
global fg = initfg()

global initCov = diagm( [1;1;1; 0.01;0.01;0.01].^2 )
global odoCov = deepcopy(initCov)

end


@testset "Adding PriorPose3 to graph..." begin

global v1 = addVariable!(fg, :x1, Pose3, N=N)
global initPosePrior = PriorPose3( MvNormal( zeros(6), initCov) )
global f1  = addFactor!(fg,[v1], initPosePrior, graphinit=true)

# ls(fg, :x1)

ensureAllInitialized!(fg)

end

@testset "Ensure variables initialized properly" begin

# start with to tight an initialization
global muX1 = mean(SpecialEuclidean(3), getVal(fg,:x1))
@test isapprox(muX1.parts[1], [0,0,0], atol=0.4)
@test isapprox(muX1.parts[2], [1 0 0; 0 1 0; 0 0 1], atol=0.04)

# global stdX1 = std(getVal(fg,:x1),dims=2)
global stdX1 = sqrt.(diag(cov(Pose3(), getVal(fg,:x1))))
@test sum(map(Int,abs.(1.0 .- stdX1[1:3]) .< 0.3)) == 3
@test sum(map(Int,abs.(0.01 .- stdX1[4:6]) .< 0.1)) == 3

end


@testset "Testing PriorPose3 evaluation..." begin

global priorpts = approxConv(fg, :x1f1, :x1)
# priorpts = evalFactor(fg, fg.g.vertices[2], 1)
global means = mean(SpecialEuclidean(3), priorpts)
@test sum(map(Int,abs.(means.parts[1]) .> 0.5)) == 0
@test isapprox(means.parts[2], [1 0 0; 0 1 0; 0 0 1], atol=0.05)

# global v2, f2 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )
# global v3, f3 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

addVariable!(fg, :x2, Pose3)
addFactor!(fg, [:x1;:x2], Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

addVariable!(fg, :x3, Pose3)
addFactor!(fg, [:x2;:x3], Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

end


@testset "Testing Pose3Pose3 evaluation..." begin

global X1pts = getVal(fg, :x1)
global X2pts = approxConv(fg, :x1x2f1, :x2, N=N)
# X2pts = evalFactor(fg, fg.g.vertices[6], 3, N=N)
global X3pts = approxConv(fg, :x2x3f1, :x3)
global X2ptsMean = mean(SpecialEuclidean(3), X2pts)
global X3ptsMean = mean(SpecialEuclidean(3), X3pts)

# @show X2ptsMean
# @show X3ptsMean

@test isapprox(X2ptsMean.parts[1], [25,0,0], atol=5.0)
@test isapprox(X2ptsMean.parts[2], [1 0 0; 0 1 0; 0 0 1], atol=0.5)

@test isapprox(X3ptsMean.parts[1], [50,0,0], atol=5.0)
@test isapprox(X3ptsMean.parts[2], [1 0 0; 0 1 0; 0 0 1], atol=0.5)

tree,smt,hist = solveTree!(fg)

end


@testset "Adding Pose3Pose3 w nullhypo to graph..." begin

global odoc3 = Pose3Pose3( MvNormal([-35;0;0;0;0;0.0], odoCov) )

@test length(intersect(DFG.listVariables(fg), [:x1;:x2;:x3])) == 3
# define 50/50% hypothesis
addFactor!(fg,[:x3;:x1],odoc3, nullhypo=0.5)

global X1pts = approxConv(fg, :x3x1f1, :x1, N=N)
global X2pts = approxConv(fg, :x3x1f1, :x3, N=N)

global p1 = manikde!(X1pts, Pose3)
global p2 = manikde!(X2pts, Pose3)

end


@testset "loading validation data for testing." begin

global tstdtdir = dirname(@__FILE__)
_X1ptst = readdlm(joinpath(tstdtdir, "X1ptst.csv"),',')
global X1ptst = map(X->getPoint(Pose3, X), eachcol(_X1ptst))
_X2ptst = readdlm(joinpath(tstdtdir, "X2ptst.csv"),',')
global X2ptst = map(X->getPoint(Pose3, X), eachcol(_X2ptst))

global p1t = manikde!(X1ptst, Pose3)
global p2t = manikde!(X2ptst, Pose3)

# plotKDE([p2t;p2],c=["red";"blue"],dims=[1;2],levels=3)
# kld(marginal(p1,[2]), marginal(p1t,[2]), method=:unscented)

global t1 = mmd(X1ptst, X1pts, AMP.SE3_Manifold, length(X1ptst), length(X1pts), bw=[0.001;])
global t2 = mmd(X2ptst, X2pts, AMP.SE3_Manifold, length(X1ptst), length(X1pts), bw=[0.001;])
# TODO Change to mmd
# global t1 = minimum([abs(kld(p1, p1t)[1]) ; abs(kld(p1t, p1)[1])])
# global t2 = minimum([abs(kld(p2, p2t)[1]) ; abs(kld(p2t, p2)[1])])

@test t1 < 0.5
@test t2 < 0.5

end



# using KernelDensityEstimatePlotting
# plotKDE([p1; p1t], dims=[1], c=["black";"red"])
