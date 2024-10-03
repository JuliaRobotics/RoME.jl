
using RoME
using Test
using DelimitedFiles
using Statistics
using TensorCast

##

# TODO this initial test is likely obsolete
@testset "Test syntax for null hypothesis matrix substitution" begin
##

function ff(X::Array{Float64,2})
  return X.+1.0
end

c = Categorical([0.5;0.5])
A = zeros(3,5)

skipdos = rand(c, 5)
dos = skipdos.==2

B = zeros(3,5)
B[:,dos] = ff(A[:,dos])

@test sum(sum(B[:,dos],dims=2) .== sum(dos),dims=1)[1,1] == 3

##
end


@testset "Test if null hypothesis occurs as expected..." begin
##

N = 100
global fg = initfg()

initCov = diagm( [1;1;1; 0.01;0.01;0.01].^2 )
odoCov = deepcopy(initCov)

# end
# @testset "Adding PriorPose3 to graph..." begin

v1 = addVariable!(fg, :x1, Pose3, N=N)
initPosePrior = PriorPose3( MvNormal( zeros(6), initCov) )
f1  = addFactor!(fg,[v1], initPosePrior, graphinit=true)

# ls(fg, :x1)

initAll!(fg)

# end
# @testset "Ensure variables initialized properly" begin

# start with to tight an initialization
muX1 = mean(SpecialEuclidean(3; vectors=HybridTangentRepresentation()), getVal(fg,:x1))
@test isapprox(submanifold_component(muX1,1), [0,0,0], atol=0.4)
@test isapprox(submanifold_component(muX1,2), [1 0 0; 0 1 0; 0 0 1], atol=0.04)

# stdX1 = std(getVal(fg,:x1),dims=2)
stdX1 = sqrt.(diag(cov(Pose3(), getVal(fg,:x1))))
@test sum(map(Int,abs.(1.0 .- stdX1[1:3]) .< 0.3)) == 3
@test sum(map(Int,abs.(0.01 .- stdX1[4:6]) .< 0.1)) == 3

# end
# @testset "Testing PriorPose3 evaluation..." begin

priorpts = approxConv(fg, :x1f1, :x1)
# priorpts = evalFactor(fg, fg.g.vertices[2], 1)
means = mean(SpecialEuclidean(3; vectors=HybridTangentRepresentation()), priorpts)
@test sum(map(Int,abs.(submanifold_component(means,1)) .> 0.5)) == 0
@test isapprox(submanifold_component(means,2), [1 0 0; 0 1 0; 0 0 1], atol=0.05)

# v2, f2 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )
# v3, f3 = addOdoFG!(fg, Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

addVariable!(fg, :x2, Pose3)
addFactor!(fg, [:x1;:x2], Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

addVariable!(fg, :x3, Pose3)
addFactor!(fg, [:x2;:x3], Pose3Pose3( MvNormal([25;0;0;0;0;0.0], odoCov)) )

##
end


@testset "Testing Pose3Pose3 evaluation..." begin
##

global fg

X1pts = getVal(fg, :x1)
X2pts = approxConv(fg, :x1x2f1, :x2, N=N)
# X2pts = evalFactor(fg, fg.g.vertices[6], 3, N=N)
X3pts = approxConv(fg, :x2x3f1, :x3)
X2ptsMean = mean(SpecialEuclidean(3; vectors=HybridTangentRepresentation()), X2pts)
X3ptsMean = mean(SpecialEuclidean(3; vectors=HybridTangentRepresentation()), X3pts)

# @show X2ptsMean
# @show X3ptsMean

@test isapprox(submanifold_component(X2ptsMean,1), [25,0,0], atol=5.0)
@test isapprox(submanifold_component(X2ptsMean,2), [1 0 0; 0 1 0; 0 0 1], atol=0.5)

@test isapprox(submanifold_component(X3ptsMean,1), [50,0,0], atol=5.0)
@test isapprox(submanifold_component(X3ptsMean,2), [1 0 0; 0 1 0; 0 0 1], atol=0.5)

tree = solveTree!(fg)

# end
# @testset "Adding Pose3Pose3 w nullhypo to graph..." begin

odoc3 = Pose3Pose3( MvNormal([-35;0;0;0;0;0.0], odoCov) )

@test length(intersect(DFG.listVariables(fg), [:x1;:x2;:x3])) == 3
# define 50/50% hypothesis
addFactor!(fg,[:x3;:x1],odoc3, nullhypo=0.5)

X1pts = approxConv(fg, :x3x1f1, :x1, N=N)
X2pts = approxConv(fg, :x3x1f1, :x3, N=N)

p1 = manikde!(Pose3, X1pts)
p2 = manikde!(Pose3, X2pts)


@info "loading validation data for testing."
tstdtdir = dirname(@__FILE__)
_X1ptst = readdlm(joinpath(tstdtdir, "X1ptst.csv"),',')
@cast _X1p[j][i] := _X1ptst[i,j]
X1ptst = map(X->DFG.getPoint(Pose3, X), _X1p)
_X2ptst = readdlm(joinpath(tstdtdir, "X2ptst.csv"),',')
@cast _X2p[j][i] := _X2ptst[i,j]
X2ptst = map(X->DFG.getPoint(Pose3, X), _X2p)

p1t = manikde!(Pose3, X1ptst)
p2t = manikde!(Pose3, X2ptst)

# plotKDE([p2t;p2],c=["red";"blue"],dims=[1;2],levels=3)

t1 = mmd(SpecialEuclidean(3; vectors=HybridTangentRepresentation()), X1ptst, X1pts, bw=[0.001;])
t2 = mmd(SpecialEuclidean(3; vectors=HybridTangentRepresentation()), X2ptst, X2pts, bw=[0.001;])

@test t1 < 0.6
@test t2 < 0.6

##
end

#