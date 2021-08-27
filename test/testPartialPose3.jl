# test partial pose3 constraints and evaluation

using Statistics
using RoME
using Test
using Manifolds: ProductRepr, hat, vee, identity_element, SpecialOrthogonal, SpecialEuclidean
import Manifolds
using TensorCast
using DistributedFactorGraphs
using Rotations

@testset "Testing basic PriorPose3ZRP" begin

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

sf = map(x->x[1], sampleFactor(fg, :x1f2, 100))
mu = getCoordinates(Pose3, mean(M, sf))

solveTree!(fg)

mpts = getPoints(fg, :x1)
mu = mean(M, mpts)
mucrd = getCoordinates(Pose3, mu)

@test isapprox(mucrd[1:3], [0, 5, 10], atol=1.0)
@test isapprox(mucrd[4:6], [0, 0, pi/2], atol=0.3)


end

# using CoordinateTransformations, Rotations

##
fg = initfg()

# @testset "Test Building a PriorPose3ZRP FG" begin
N = 50
fg.solverParams.N = N

v1 = addVariable!(fg,:x1, Pose3) # 0.001*randn(6,N)
# f0 = addFactor!(fg, [:x1;], PriorPose3(MvNormal(zeros(6),1e-2*Matrix{Float64}(LinearAlgebra.I, 6,6))))
# check with a point not a identity
f0 = addFactor!(fg, [:x1;], PriorPose3(MvNormal([0.0, 5.0, 10, pi, 0, 0],1e-2*Matrix{Float64}(LinearAlgebra.I, 6,6))))

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

@testset "test PriorPose3ZRP evaluations" begin

##

# ensure that at least the first pose is already initialized
doautoinit!(fg, :x1)
@test isInitialized(fg, :x1)

X1pts = getVal(fg, :x1)
# @test sum(isnan.(X1pts)) == 0 
# ProductRepr does not have a NaN, so testing for ProductRepr instead
@test all(isa.(X1pts, ProductRepr))

ppts = approxConv(fg, :x1x2f1, :x2)
@test all(isa.(ppts, ProductRepr))


initAll!(fg)
@test isInitialized(fg, :x2)

# get values and ensure that a re-evaluation produces consistent results
_X2pts = getCoordinates.(Pose3, getVal(fg, :x2))
@cast X2pts[j,i] :=  _X2pts[i][j]
@test sum(isnan.(X2pts)) == 0

# check that only partial states are updated
_pts = getCoordinates.(Pose3, IIF.approxConv(fg, :x2f1, :x2, N=N))

@cast pts[j,i] :=  _pts[i][j]

newdims = collect(DFG.getSolverData(f1).fnc.usrfnc!.partial)

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

tfg = initfg()
X0 = addVariable!(tfg, :x0, Pose3)
X1 = addVariable!(tfg, :x1, Pose3)

M = getManifold(Pose3)
ϵ = identity_element(M)
xi = deepcopy(ϵ)
xja = deepcopy(ϵ)

res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), Tuple[], (xi, xja))


@test abs(res[1]-mu2[1]) < 0.3
@test abs(res[2]-mu2[2]) < 0.3
@test abs(res[3]-mu2[3]) < 0.2

##

Xjb = zeros(6,1)
Xjb[collect(xyy.partial),1] = mu2
xjb = exp(M, ϵ, hat(M, ϵ, Xjb))
# res = zeros(3)
# xyy(res, fmd, idx, meas, xi, xjb)

res = calcFactorResidualTemporary(xyy, (Pose3, Pose3), Tuple[], (xi, xjb))

@test 0.0 < norm(res) < 0.3

##

addFactor!(tfg, [:x0;:x1], xyy, graphinit=false)

# meas = getSample(xyy,100)
ccw = IIF._getCCW(tfg, :x0x1f1)
_meas = sampleFactor(ccw, 100)

meas = map(x->x[1], _meas)

M = getManifold(xyy)

@test isapprox(sqrt.(diag(cov(M, meas))), [0.01;0.01;0.002], atol=0.05)

##

end



@testset "test Pose3Pose3XYYaw evaluations" begin

##

# get existing and predict new
_X2pts = getCoordinates.(Pose3, getBelief(fg, :x2) |> getPoints)
@cast X2pts[j,i] :=  _X2pts[i][j]

_pts = getCoordinates.(Pose3, approxConv(fg, :x1x2f1, :x2, N=N))
@cast pts[j,i] :=  _pts[i][j]

pts = collect(pts)

# find which dimensions are and and are not updated by XYYaw partial
newdims = collect(DFG.getSolverData(f2).fnc.usrfnc!.partial)
olddims = setdiff(collect(1:6), newdims)

# check the number of points are correct
@test size(pts, 1) == 6
@test size(pts, 2) == N

# ensure the unchanged dimensions actually remain unchanged
@test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

for i in 1:N
    pts[6,i] = wrapRad(pts[6,i])
end

# ensure the newly updated values match what is specified in mu2
@show mu2 # mu2 is used for XYYaw
# trying to compare world and body frame values -- not right!
@show Statistics.mean(pts[newdims,:],dims=2)
@test sum(abs.(Statistics.mean(pts[newdims,:],dims=2)-mu2) .< [1.5;1.5;0.3]) == 3

# ensure a re-evaluation of the partial factor updates the partial variable dimensions correclty
@test norm(X2pts[newdims,:] - pts[newdims,:]) < 1.0

# ensure that memory pointers are working correctly
memcheck = getVal(v2)
# @test norm(X2pts - memcheck) < 1e-10

##

end


@testset "test predictbelief with two functions" begin

##

_val = getCoordinates.(Pose3, predictbelief(fg, :x2, ls(fg, :x2), N=N)[1])
@cast val[j,i] :=  _val[i][j]
val = copy(val)
for i in 1:N
  val[6,i] = wrapRad(val[6,i])
end

@test size(val, 1) == 6
@test size(val, 2) == N

estmu1mean = Statistics.mean(val[collect(DFG.getSolverData(f1).fnc.usrfnc!.partial),:],dims=2)
estmu2mean = Statistics.mean(val[collect(DFG.getSolverData(f2).fnc.usrfnc!.partial),:],dims=2)

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
  [ 0.,  0, 10, 0, 0, ψ],
]

@info "Test Pose3Pose3XYYaw cases:"
xyz_rpy = testsMeasurements_xyz_rpy[1]
for xyz_rpy = testsMeasurements_xyz_rpy

@info xyz_rpy
wTx1 = getPointIdentity(Pose3)
wTx2 = ProductRepr(xyz_rpy[1:3], RotXYZ(xyz_rpy[4:6]...))

x1Tx2 = Manifolds.compose(M3, inv(M3, wTx1) , wTx2)
X = log(M3, ϵ3, x1Tx2)
Xc = vee(M3, ϵ3, X)

XYH1_2 = [xyz_rpy[1], xyz_rpy[2], xyz_rpy[6]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(Xc, 0.001*diagm(ones(6))))

meas = [(X,)]

res = calcFactorResidualTemporary(testpp3, (Pose3, Pose3), meas, (wTx1, wTx2))

@test norm(res) < 1e-10


# test with PartialXYH
testppxyh = Pose3Pose3XYYaw( MvNormal(XYH1_2, 0.001*diagm([1.0;1;1])))

Xc = hat(M2, ϵ2, XYH1_2)
meas = [(Xc, )]
res = calcFactorResidualTemporary(testppxyh, (Pose3, Pose3), meas, (wTx1, wTx2))


@test norm(res) < 1e-10

end
##

end

