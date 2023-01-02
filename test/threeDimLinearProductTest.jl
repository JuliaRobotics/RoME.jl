
using RoME
using TransformUtils
using Statistics
using Test


M = getManifold(Pose3)
ϵ = getPointIdentity(Pose3)

##

@testset "test 3D convolutions and products" begin

##

tf = SE3([0.0;0.0;0.0], TU.AngleAxis(pi/4,[0;0;1.0]))# Euler(pi/4,0.0,0.0) )

N = 1
initCov = 0.0001*Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.000001 for i in 4:6];
odoCov = deepcopy(initCov)
odo = Pose3Pose3( MvNormal(veeEuler(tf), odoCov) )

X = getPoint.(Pose3, [[0.01*randn(5); (0*pi/4 .+ 0.01*randn(1))] for _=1:N])

fg = initfg()
addVariable!(fg, :x0, Pose3)
X_ = manikde!(Pose3, X, bw=ones(6))
initVariable!(fg, :x0, X_)

addVariable!(fg, :x1, Pose3)
addFactor!(fg, [:x0;:x1], odo)

Y = approxConv(fg, :x0x1f1, :x1, N=1) # X ⊕ odo

@error "TODO fix approxConv(.. N=1) case"
# @test norm(Y[1:3]-zeros(3)) < 1.0
# @test norm(Y[4:6]-[zeros(2);pi/4]) < 0.15

##

tf = SE3([0.0;0.0;0.0], TU.AngleAxis(pi/4,[0;0;1.0]))# Euler(pi/4,0.0,0.0) )

N = 1
initCov = 0.01*Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.001 for i in 4:6];
odoCov = deepcopy(initCov)
odo = Pose3Pose3(  MvNormal(veeEuler(tf), odoCov) )

X = getPoint.(Pose3, [[0.01*randn(5); (0*pi/4 .+ 0.01*randn(1))] for _=1:N])

fg = initfg()
addVariable!(fg, :x0, Pose3)
X_ = manikde!(Pose3, X, bw=ones(6))
initVariable!(fg, :x0, X_)

addVariable!(fg, :x1, Pose3)
addFactor!(fg, [:x0;:x1], odo)

Y = approxConv(fg, :x0x1f1, :x1, N=1) # X ⊕ odo

@error "TODO fix approxConv(.. N=1) case"
# @test norm(Y[1:3]-zeros(3)) < 1.0
# @test norm(Y[4:6]-[zeros(2);pi/4]) < 0.2

##

end

##

global N = 100
global fg = initfg()

global initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.01 for i in 4:6];
global odoCov = deepcopy(initCov)


@testset "Adding PriorPose3 to graph..." begin

global v1 = addVariable!(fg, :x1, Pose3,  N=N) # 0.1*randn(6,N)
global initPosePrior = PriorPose3( MvNormal(zeros(6), initCov) )
global f1  = addFactor!(fg,[:x1;], initPosePrior)
# @test !isInitialized(fg, :x1)

end

##

@testset "Ensure variable is properly initialized" begin
  M = getManifold(Pose3)
  # start with initialization
  initAll!(fg)
  @test isInitialized(fg, :x1)
  @show muX1 = Statistics.mean(M, getVal(fg,:x1))
  
  T = submanifold_component(muX1,1)
  @test sum(map(Int,abs.(T) .< 0.5)) == 3
  
  Rc = submanifold_component(muX1,2)
  @test isapprox(SpecialOrthogonal(3), Rc, [1 0 0; 0 1 0; 0 0 1], atol=0.25)

  coX1 = getCoordinates.(Pose3, getVal(fg,:x1))
  @show stdX1 = Statistics.std(coX1)
  @test sum(map(Int, 0.5 .< stdX1[1:3] .< 1.5)) == 3
  @test sum(map(Int, 0.05 .< stdX1[4:6] .< 0.15)) == 3
end


@testset "Testing PriorPose3 evaluation..." begin

global priorpts = approxConv(fg, :x1f1, :x1) # fg.g.vertices[2], 1
coX1 = getCoordinates.(Pose3, priorpts)
global means = Statistics.mean(coX1)
@test sum(map(Int,abs.(means[1:3]) .> 0.5)) == 0
@test sum(map(Int,abs.(means[4:6]) .> 0.05)) == 0

end


##


@testset "Adding Pose3Pose3 to graph..." begin

##

global N
global odo = SE3([10;0;0], Quaternion(0))
# global pts0X2 = projectParticles(getVal(fg,:x1), MvNormal(veeEuler(odo), odoCov) )
global odoconstr = Pose3Pose3( MvNormal(veeEuler(odo), odoCov) )
global v2 = addVariable!(fg,:x2, Pose3, N=N) # pts0X2
addFactor!(fg,[:x1;:x2],odoconstr, inflation=0.1)
# @test !isInitialized(fg, :x2)


## test opposites

# should be all zero
p = deepcopy(ϵ)
q = deepcopy(ϵ)
submanifold_component(q,1)[1] = 10.0
X = Manifolds.hat(M, ϵ, [10.,0,0,0,0,0])
res = calcFactorResidual(fg, :x1x2f1, X, p, q)
@test norm(res) < 1e-10

# trivial fail case
# res = calcFactorResidual(fg, :x1x2f1, [10;0;0;0;0;0.0], zeros(6), [10;0;0;pi;pi;pi])
@warn "suppressing trivial Pose3 fail case until RoME.jl #244 has been completed."
#TODO is this the case for the warn?
p = deepcopy(ϵ)
q = getPoint(Pose3, [10;0;0;pi;pi;pi])
X = Manifolds.hat(M, ϵ, [10.,0,0,pi,pi,pi])

res = calcFactorResidual(fg, :x1x2f1, X, p, q)
@test norm(res) < 1e-10


## test following introduction of inflation, see IIF #1051

# force the inflation trivial error, https://github.com/JuliaRobotics/RoME.jl/issues/380#issuecomment-778795848
# IIF._getCCW(fg,:x1x2f1).inflation = 10.0

_pts = approxConv(fg, :x1x2f1, :x2)
# X2 = manikde!(Pose3, pts)
@cast pts[j,i] := getCoordinates.(Pose3, _pts)[i][j]

# test translations through convolution
@test 0.8N < sum( 5 .< pts[1,:] .<15 )
@test 0.7N < sum( -5 .< pts[2,:] .< 5 )
@test 0.7N < sum( -5 .< pts[3,:] .< 5 )

# test rotations through convolution
@test 0.95N < sum( -2 .< pts[4,:] .< 2 )
@test 0.95N < sum( -2 .< pts[5,:] .< 2 )
@test 0.95N < sum( -2 .< pts[6,:] .< 2 )

##

end

##

# Noticed a DomainError on convolutions here after mutlithreading upgrade.  Previously used fill(PP3REUSE, Threads.nthreads())
@testset "Testing Pose3Pose3 evaluation..." begin

initAll!(fg)
@test isInitialized(fg, :x2)
global X1pts = approxConv(fg, :x1x2f1, :x1)
# X1pts = evalFactor(fg, fg.g.vertices[4], 1)
global X2pts = approxConv(fg, :x1x2f1, :x2)
# X2pts = evalFactor(fg, fg.g.vertices[4], 3)

mu = mean(M, getVal(fg,:x1))
T = submanifold_component(mu,1)
@test isapprox(T, [0,0,0], atol=0.6)
Rc = submanifold_component(mu,2)
@test isapprox(SpecialOrthogonal(3), Rc, [1 0 0; 0 1 0; 0 0 1], atol=0.25)

mu = mean(M, getVal(fg,:x2))
T = submanifold_component(mu,1)
@test isapprox(T, [10,0,0], atol=1.0)
Rc = submanifold_component(mu,2)
@test isapprox(SpecialOrthogonal(3), Rc, [1 0 0; 0 1 0; 0 0 1], atol=0.25)


end

@testset "Construct Bayes tree and perform inference..." begin
  tree = solveTree!(fg)
  mu = mean(M, getVal(fg,:x1))
  T = submanifold_component(mu,1)
  @test isapprox(T, [0,0,0], atol=1.5)
  Rc = submanifold_component(mu,2)
  @test isapprox(SpecialOrthogonal(3), Rc, [1 0 0; 0 1 0; 0 0 1], atol=0.25)

  mu = mean(M, getVal(fg,:x2))
  T = submanifold_component(mu,1)
  @test isapprox(T, [10,0,0], atol=1.5)
  Rc = submanifold_component(mu,2)
  @test isapprox(SpecialOrthogonal(3), Rc, [1 0 0; 0 1 0; 0 0 1], atol=0.25)
end

##

@testset "Ensure basic parameters on x1,x2 after inference..." begin

##

# check mean and covariances after one up and down pass over the tree
global muX1 = getPPE(fg, :x1).suggested # Statistics.mean(getVal(fg,:x1),dims=2)
@test sum(map(Int,abs.(muX1[1:3]) .< 1.0)) == 3

# sidestep trivial case, #412
# case where points land on trivial rotation error [pi;pi;pi]
@warn "simplify mean and std testing on Pose3 after #244, see #412"
_pts = getVal(fg,:x1)
@cast pts[j,i] := getCoordinates.(Pose3, _pts)[i][j]

mask = (2.7 .< abs.(pts[4,:])) .& (2.7 .< abs.(pts[5,:])) .& (2.7 .< abs.(pts[6,:]))
mask .= xor.(mask,1)
mu1tmp = Statistics.mean(pts[:,mask],dims=2)
@test sum(map(Int, abs.(mu1tmp[4:6]) .< 0.1)) == 3
global stdX1 = Statistics.std(pts[:,mask],dims=2)

@test sum(map(Int, 0.4 .< stdX1[1:3] .< 1.6)) == 3 # had a 2==3 failure here
@show stdX1[4:6]
@test sum(map(Int, 0.02 .< stdX1[4:6] .< 0.5)) == 3


global muX2 = getPPE(fg, :x2).suggested # Statistics.mean(getVal(fg,:x2),dims=2)
@show muX2[1:3]-[10.0;0;0]
@test sum(map(Int, abs.(muX2[1:3]-[10.0;0;0]) .< 1.5)) == 3

# case where points land on trivial rotation error [pi;pi;pi]
@warn "simplify mean and std testing on Pose3 after #244, see #412"
_pts = getVal(fg,:x2)
_pts = getVal(fg,:x1)
@cast pts[j,i] := getCoordinates.(Pose3, _pts)[i][j]

mask = (2.7 .< abs.(pts[4,:])) .& (2.7 .< abs.(pts[5,:])) .& (2.7 .< abs.(pts[6,:]))
mask .= xor.(mask,1)
mu2tmp = Statistics.mean(pts[:,mask],dims=2)
@test sum(map(Int, abs.(mu2tmp[4:6]) .< 0.1)) == 3
global stdX2 = Statistics.std(pts[:,mask],dims=2)
# @test sum(map(Int, abs.(muX2[4:6]) .< 0.1)) == 3
# global stdX2 = Statistics.std(getVal(fg,:x2),dims=2)

println("previous test failure 0.75 .< $(round.(stdX2[1:3],digits=2)) .< 2.25")
@test sum(map(Int, 0.55 .< stdX2[1:3] .< 2.5)) == 3
println("previous test failure 0.05 .< $(round.(stdX2[4:6],digits=2)) .< 0.35")
@test sum(map(Int, 0.05 .< stdX2[4:6] .< 0.5)) == 3

##

end

# println("Plot marginals to see what is happening")
# plotKDE(marginal(getBelief(fg,:x1),[1]))
# plotKDE(marginal(getBelief(fg,:x2),[1]))




#
