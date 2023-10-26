# test partial pose3 constraints and evaluation

# using Revise
using Statistics
using RoME
using Test
using Manifolds: hat, vee, identity_element, SpecialOrthogonal, SpecialEuclidean
import Manifolds
using TensorCast
using DistributedFactorGraphs
using Rotations

##


@testset "Test (partial) priors on Pose3" begin
## there was a problem with partials and staticarrays where asigning PriorPoint3 to Pose3 stopped working, here is an explicit test

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :x0, Pose3)

mv3 = MvNormal(SA[0.0; 0.0; 0.0], SMatrix{3,3}(1,0,0,0,1,0,0,0,1.))
f = addFactor!(fg, [:x0;], PriorPoint3(mv3))

approxConvBelief(fg, getLabel(f), :x0)

##
end


@testset "Testing basic PriorPose3ZRP" begin
##

fg = initfg()

M=SpecialEuclidean(3)
N = 100
fg.solverParams.N = N
fg.solverParams.graphinit = false

v1 = addVariable!(fg,:x1, Pose3)
#                                                 x    y    z    ϕ    θ   ψ
f0 = addFactor!(fg, [:x1], PriorPose3(MvNormal([0.0, 5.0, 9.0, 0.1, 0.0, pi/2], diagm([1,1,1,0.1,0.1,0.1].^2) )))
#                               z                      ϕ    θ  
prpz = PriorPose3ZRP( Normal(11.0, 1.0), MvNormal( [-0.1, 0.0], diagm([0.1, 0.1].^2) ))
f1 = addFactor!(fg, [:x1], prpz)

sf = sampleFactor(fg, :x1f2, N)

Mzrp = f1 |> getFactorType |> getManifold
zrp0 = identity_element(Mzrp)
mu = vee(Mzrp, zrp0, log(Mzrp, zrp0, mean(Mzrp, sf)))
# mu = getCoordinates(Pose3, mean(Mzrp, sf)) # M # starting to improve partials, getManifold of partial factor returns whatever that manifold is

solveTree!(fg)

mpts = getPoints(fg, :x1)
mu = mean(M, mpts)
mucrd = getCoordinates(Pose3, mu)

@test isapprox(mucrd[1:3], [0, 5, 10], atol=1.5)
@test isapprox(mucrd[4:6], [0, 0, pi/2], atol=0.3)

##
end

##

fg = initfg()

# @testset "Test Building a PriorPose3ZRP FG" begin
N = 50
fg.solverParams.N = N

v1 = addVariable!(fg,:x1, Pose3) # 0.001*randn(6,N)
# f0 = addFactor!(fg, [:x1;], PriorPose3(MvNormal(zeros(6),1e-2*Matrix{Float64}(LinearAlgebra.I, 6,6))))
# check with a point not a identity
f0 = addFactor!(fg, [:x1;], PriorPose3(MvNormal([0.0, 0.0, 10, 0, 0, 0],1e-2*diagm(ones(6)))))

sigx = 0.01
sigy = 0.01
sigth = 0.0281
mu1 = [0.0;0.0; -10.0]
prpz = PriorPose3ZRP(
  Normal( mu1[3], sigth ),
  MvNormal( mu1[1:2], [sigx 0.0; 0.0 sigy].^2 )
)

mu2 = [20.0,5.0,pi/2]
xyy = Pose3Pose3XYYaw(MvNormal( mu2, diagm([sigx, sigy, sigth].^2)))

v2 = addVariable!(fg,:x2, Pose3)

f1 = addFactor!(fg, [:x2], prpz, graphinit=false)
f2 = addFactor!(fg, [:x1;:x2], xyy, graphinit=false)

#
# end

##

@testset "test PriorPose3ZRP evaluations" begin
##

# ensure that at least the first pose is already initialized
doautoinit!(fg, :x1)
@test isInitialized(fg, :x1)

X1pts = getVal(fg, :x1)
# @test sum(isnan.(X1pts)) == 0 
# ProductRepr does not have a NaN, so testing for ProductRepr instead
@test all(isa.(X1pts, ArrayPartition))

ppts = approxConv(fg, :x1x2f1, :x2)
@test all(isa.(ppts, ArrayPartition))


initAll!(fg)
@test isInitialized(fg, :x2)

# get values and ensure that a re-evaluation produces consistent results
_X2pts = getCoordinates.(Pose3, getVal(fg, :x2))
@cast X2pts[j,i] :=  _X2pts[i][j]
@test sum(isnan.(X2pts)) == 0

# check that only partial states are updated
_pts = getCoordinates.(Pose3, IIF.approxConv(fg, :x2f1, :x2, N=N))

@cast pts[j,i] :=  _pts[i][j]

newdims = [3;4;5] # collect(DFG.getSolverData(f1).fnc.usrfnc!.partial)

olddims = setdiff(collect(1:6), newdims)


@test size(pts, 1) == 6
@test size(pts, 2) == N

# check that untouched dimensions (not in the partial list) truely remain untouched
@test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

# check that the prior new dims are updated to new and correct values
# @show Statistics.mean(pts,dims=2)[newdims]
@test sum(abs.(Statistics.mean(pts,dims=2)[newdims]-mu1[[3;1;2]]) .< [0.5; 0.1; 0.1]) == 3

# ensure a forced re-evaluatoin
@test norm(X2pts[newdims,:] - pts[newdims,:]) < 1.0

# memcheck that the exact same values are used
# @test norm(X2pts - getVal(v2)) < 1e-10

##
end


@testset "test residual function of Pose3Pose3XYYaw" begin
##

M = getManifold(Pose3)
ϵ = getPointIdentity(M)
xi = deepcopy(ϵ)
xja = deepcopy(ϵ)

res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xja))

@test isapprox(mu2, res, atol=0.3 )

##

Xjb = zeros(6,1)
Xjb[[1;2;6],1] = mu2  
# Xjb[collect(xyy.partial),1] = mu2

# noisy measurements
xi = ArrayPartition([0;0;0.], collect(RotZ(0)))
xjb = exp(M, ϵ, hat(M, ϵ, Xjb))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xjb))
@test isapprox(res, [0;0;0], atol=0.15)

# more rotated
xi = ArrayPartition([0;0;0.], collect(RotZ(π/2)))
xj = ArrayPartition([-5;20;0.], [-1 0 0; 0 -1 0; 0 0 1.])
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add z
xi = ArrayPartition([0;0;100.], collect(RotZ(π/2)))
xj = ArrayPartition([-5;20;-100.], [-1 0 0; 0 -1 0; 0 0 1.])
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add pitch without z
xi = ArrayPartition([0;0;0.], collect(RotY(π/4)))
xj = ArrayPartition([20;5;0.], collect(RotZ(π/2)*RotY(π/4)))
# xj = ArrayPartition([14.14213562;5;-14.14213562], [0 -0.707107 0.707107; 1 0 0; 0 0.707107 0.707107] )
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add pitch with z
xi = ArrayPartition([0;0;10.], collect(RotY(π/4)))
xj = ArrayPartition([20;5;-10.], collect(RotZ(π/2)*RotY(π/4)))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add roll without z
xi = ArrayPartition([0;0;0.], collect(RotX(π/4)))
xj = ArrayPartition([20;5;0.], collect(RotZ(π/2)*RotX(π/4)))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add roll with z
xi = ArrayPartition([0;0;10.], collect(RotX(π/4)))
xj = ArrayPartition([20;5;-10.], collect(RotZ(π/2)*RotX(π/4)))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add roll with x
xi = ArrayPartition([10;0;0.], collect(RotX(π/4)))
xj = ArrayPartition([10+20;5;0.], collect(RotZ(π/2)*RotX(π/4)))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add pitch and roll without z
xi = ArrayPartition([0;0;0.], collect(RotY(π/4)*RotX(π/6)))
xj = ArrayPartition([20;5;0.], collect(RotZ(π/2)*RotY(π/4)*RotX(π/6)))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)

# add pitch and roll with x and z
xi = ArrayPartition([10;0;10.], collect(RotY(π/4)*RotX(π/6)))
xj = ArrayPartition([10+20;5;-10.], collect(RotZ(π/2)*RotY(π/4)*RotX(π/6)))
res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [0;0;0], atol=0.15)



##

tfg = initfg()
X0 = addVariable!(tfg, :x0, Pose3)
X1 = addVariable!(tfg, :x1, Pose3)

addFactor!(tfg, [:x0;:x1], xyy, graphinit=false)

# meas = getSample(xyy,100)
ccw = IIF._getCCW(tfg, :x0x1f1);
meas = sampleFactor(ccw, 100)

M = getManifold(xyy)

@test isapprox(sqrt.(diag(cov(M, meas))), [0.01;0.01;0.002], atol=0.05)

## Testing non zero residuals for sign
# Lets just make sure we are consistent with signs
lr = LinearRelative(Normal(20.0,0.01))
xi = [10]
xj = [20]
res = calcFactorResidualTemporary(lr, (ContinuousScalar, ContinuousScalar), [], (xi, xj))
@test isapprox(res, [10], atol=0.1)

# To compare signs to Point2 case 
mu = [20.0, 5.0]
xyy2 = Point2Point2(MvNormal( mu, diagm([0.01, 0.01].^2)))
xi = [10;0]
xj = [20;10]
res = calcFactorResidualTemporary(xyy2, (Point2, Point2), [], (xi, xj))
@test isapprox(res, [10;-5], atol=0.15)

# Comparing to Point2 for sign only
mu = [20.0, 5.0, 0.0]
xyy2 = Pose3Pose3XYYaw(MvNormal( mu, diagm([0.01, 0.01, 0.001].^2)))
xi = ArrayPartition([10;0; 10.], collect(RotZ(0.0)))
xj = ArrayPartition([20;10;-10.], collect(RotZYX(0.0, -pi/4, pi/4)))
res = calcFactorResidualTemporary(xyy2, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [10;-5;0], atol=0.15)


mu = [20.0, 5.0, pi/4]
xyy2 = Pose3Pose3XYYaw(MvNormal( mu, diagm([0.01, 0.01, 0.001].^2)))
xi = ArrayPartition([10;0; 10.], collect(RotZ(pi/2)))
xj = ArrayPartition([20;10;-10.], collect(RotZYX(0.0, -pi/4, pi/4)))
res = calcFactorResidualTemporary(xyy2, (Pose3, Pose3), [], (xi, xj))
@test isapprox(res, [-15;10;3pi/4], atol=0.15)

##
end



@testset "test Pose3Pose3XYYaw evaluations" begin
##

# get existing and predict new
_X2 = getBelief(fg, :x2) 
_X2pts_ = _X2 |> getPoints

_X2prd = approxConvBelief(fg, :x1x2f1, :x2, N=N)

# test for #1394
@test isPartial(_X2prd)

_X2prd_ = getPoints(_X2prd, false)

@test isapprox(submanifold_component(mean(_X2prd, false),1), submanifold_component(mean(_X2, false),1), atol=2.0)
@test isapprox(submanifold_component(mean(_X2prd, false),2), submanifold_component(mean(_X2, false),2), atol=0.5)

convert(TU.Euler, SO3(collect(submanifold_component(mean(_X2prd, false),2)))).Y
convert(TU.Euler, SO3(collect(submanifold_component(mean(_X2prd, false),2)))).P
convert(TU.Euler, SO3(collect(submanifold_component(mean(_X2prd, false),2)))).R

convert(TU.Euler, SO3(collect(submanifold_component(mean(_X2, false),2)))).Y
convert(TU.Euler, SO3(collect(submanifold_component(mean(_X2, false),2)))).P
convert(TU.Euler, SO3(collect(submanifold_component(mean(_X2, false),2)))).R


# pts = collect(pts)

# find which dimensions are and and are not updated by XYYaw partial
newdims = [1,2,6] # collect(DFG.getSolverData(f2).fnc.usrfnc!.partial)
olddims = setdiff(collect(1:6), newdims)

# check the number of points are correct
# @test size(pts, 1) == 6
@test length(_X2pts_) == N


## DEV

i = 1
for i in 1:N

@test isapprox( submanifold_component(_X2prd_[i],1)[3], submanifold_component(_X2pts_[i],1)[3], atol=0.0001 )

@test isapprox( convert(TU.Euler, SO3(collect(submanifold_component(_X2prd_[i],2)))).R,
                convert(TU.Euler, SO3(collect(submanifold_component(_X2pts_[i],2)))).R, atol=0.2 )

@test isapprox( convert(TU.Euler, SO3(collect(submanifold_component(_X2prd_[i],2)))).P,
                convert(TU.Euler, SO3(collect(submanifold_component(_X2pts_[i],2)))).P, atol=0.2 )

end

# ensure the unchanged dimensions actually remain unchanged

##
end


@testset "test predictbelief with two functions" begin
##

_val = IIF.getCoordinates.(Pose3, getPoints(propagateBelief(fg, :x2, ls(fg, :x2), N=N)[1]))
@cast val[j,i] :=  _val[i][j]
val = copy(val)
for i in 1:N
  val[6,i] = wrapRad(val[6,i])
end

@test size(val, 1) == 6
@test size(val, 2) == N

estmu1mean = Statistics.mean(val[collect(DFG.getSolverData(f1).fnc.usrfnc!.partial),:],dims=2)
# estmu2mean = Statistics.mean(val[collect(DFG.getSolverData(f2).fnc.usrfnc!.partial),:],dims=2)
estmu2mean = Statistics.mean(val[[1,2,6],:],dims=2)

@show estmu1mean
@show estmu2mean
@test sum(abs.(estmu1mean - mu1[[3;1;2]]) .< [0.7; 0.1; 0.1]) == 3
@test sum(abs.(estmu2mean - mu2) .< [1.0; 1.5; 0.2] ) == 3

memcheck = getVal(v2)
# @test 1e-10 < norm(val - memcheck)

##
end


# from testTartialPose3XYYaw


@testset "Test Pose3Pose3XYYaw combinations:" begin

##
# trivial cases first, orientation based tests below
M3 = getManifold(Pose3)
ϵ3 = getPointIdentity(Pose3)

M2 = getManifold(Pose2)
ϵ2 = getPointIdentity(Pose2)

ϕ = 0.1
θ = 0.1
ψ = 0.1

#
#NOTE Test is built in a way that only one angle can be non zero at a time
testsMeasurements_xyz_rpy = [
  [10.,  0,  0, 0, 0, 0],
  [ 0., 10,  0, 0, 0, 0],
  [ 0.,  0, 10, 0, 0, 0],
  [10.,  0,  0, ϕ, 0, 0],
  [ 0., 10,  0, ϕ, 0, 0],
  [ 0.,  0, 10, ϕ, 0, 0],
  [10.,  0,  0, 0, θ, 0],
  [ 0., 10,  0, 0, θ, 0],
  [ 0.,  0, 10, 0, θ, 0],
  [ 0., 15, 10, 0, θ, 0],
  [10.,  0,  0, 0, 0, ψ],
  [ 0., 10,  0, 0, 0, ψ],
  [ 0.,  0, 10, 0, 0, ψ],
  [ 0., 15, 10, 0, 0, ψ],
]

@info "Test Pose3Pose3XYYaw cases:"
xyz_rpy = testsMeasurements_xyz_rpy[1]
for xyz_rpy = testsMeasurements_xyz_rpy

@info xyz_rpy
wTx1 = getPointIdentity(Pose3)
wTx2 = ArrayPartition(xyz_rpy[1:3], Matrix(RotXYZ(xyz_rpy[4:6]...)))

x1Tx2 = Manifolds.compose(M3, inv(M3, wTx1) , wTx2)
X = log(M3, ϵ3, x1Tx2)
Xc = vee(M3, ϵ3, X)

XYH1_2 = [xyz_rpy[1], xyz_rpy[2], xyz_rpy[6]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(Xc, 0.001*diagm(ones(6))))

meas = X

res = calcFactorResidualTemporary(testpp3, (Pose3, Pose3), meas, (wTx1, wTx2))

@test norm(res) < 1e-10


# test with PartialXYH
testppxyh = Pose3Pose3XYYaw( MvNormal(XYH1_2, 0.001*diagm([1.0;1;1])))

Xc = hat(M2, ϵ2, XYH1_2)
meas = Xc
res = calcFactorResidualTemporary(testppxyh, (Pose3, Pose3), meas, (wTx1, wTx2))


@test norm(res) < 1e-10

end
##

end


@testset "Square PriorPose3ZRP and Pose3Pose3XYYaw combination solve" begin
##
# Another check

## Build 2 test fg one with σRP = 0.0 and one with σRP = 1.0

fg = (initfg(), initfg())
σRPs = [0.0 1.0]

for j=1:2
  σRP = σRPs[j]
  
  N = 4

  for i = 0:N
    addVariable!(fg[j],Symbol("x$i"), Pose3)
  end
  #Pprior                                           x    y    z    ϕ    θ    ψ
  f0 = addFactor!(fg[j], [:x0], PriorPose3(MvNormal([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], diagm([0.1,0.1,0.1,0.01,0.01,0.01].^2) )))

  for i = 1:N
    prpz = PriorPose3ZRP( Normal(i, 0.1), MvNormal( σRP*randn(2), diagm([0.01, 0.01].^2) ))
    addFactor!(fg[j], [Symbol("x$i")], prpz)
  end

  for i = 1:N
    xyy = Pose3Pose3XYYaw(MvNormal( [10.0, 0, pi/2], diagm([0.1, 0.1, 0.01].^2)))
    addFactor!(fg[j], [Symbol("x$(i-1)"), Symbol("x$i")], xyy)
  end

end


##
M = SpecialEuclidean(3)
mpts = getPoints(fg[1], :x4)
mu_fg1 = mean(M, mpts)

@test isapprox(submanifold_component(mu_fg1,1), [0,0,4], atol=0.3)

mpts = getPoints(fg[2], :x4)
mu_fg2 = mean(M, mpts)

@test_broken isapprox(submanifold_component(mu_fg2,1), [0,0,4], atol=0.3)

##
end


@testset "Test Pose3Pose3Rotation factor:" begin

##
# trivial cases first, orientation based tests below
M3 = getManifold(Pose3)
ϵ3 = getPointIdentity(Pose3)

M2 = SpecialOrthogonal(3)
ϵ2 = getPointIdentity(M2)

ϕ = 0.1
θ = 0.1
ψ = 0.1

#
#NOTE Test is built in a way that only one angle can be non zero at a time
testsMeasurements_rpy = [
  [ 0,  0,  0],
  [ ϕ,  0,  0],
  [ 0,  θ,  0],
  [ 0,  0,  ψ],
  [-ϕ,  0,  0],
  [ 0, -θ,  0],
  [ 0,  0, -ψ],
]

@info "Test Pose3Pose3Rotation cases:"
rpy = testsMeasurements_rpy[1]
for rpy = testsMeasurements_rpy

@info rpy
wTx1 = getPointIdentity(Pose3)
wTx2 = ArrayPartition([0.,0,0], Matrix(RotXYZ(rpy...)))

# test with Pose3Pose3Rotation
testppxyh = Pose3Pose3Rotation( MvNormal(rpy, 0.001*diagm([1;1;1])))

meas = hat(M2, ϵ2, rpy)
res = calcFactorResidualTemporary(testppxyh, (Pose3, Pose3), meas, (wTx1, wTx2))


@test norm(res) < 1e-10

end
##

end
