# test partial pose3 constraints and evaluation

using Statistics
using RoME
using Test
using Manifolds: ProductRepr
import Manifolds
using TensorCast
using DistributedFactorGraphs

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


## from testTartialPose3XYYaw

#TODO Update the next tests

@test_skip begin

@testset "test x translation case" begin

##
#TODO update to use manifolds and remove TransformUtils
# trivial cases first, orientation based tests below
wTx = Vector{AffineMap}(undef, 2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(UnitQuaternion(1.0, 0, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(UnitQuaternion(1.0, 0, 0, 0))

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*diagm(ones(6))))
# res = zeros(6)


# dummy value
M = Manifolds.SpecialEuclidean(3)
ϵ = identity_element(M)
meas = (Manifolds.hat(M, ϵ, veeEuler(x1Tx2)),)
p = exp(M, ϵ, Manifolds.hat(M, ϵ, veeEuler(wTx1)))
q = exp(M, ϵ, Manifolds.hat(M, ϵ, veeEuler(wTx2)))
res = calcFactorResidualTemporary(testpp3, (Pose3, Pose3), meas, (p, q))

@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

#

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(
  MvNormal(XYH1_2[1:2], 0.001*diagm([1.0;1.0])),
  Normal(XYH1_2[3], 0.001)
)
# testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2,0.001*Matrix{Float64}(LinearAlgebra.I, 3,3))  )

# res = zeros(3)
# cfo = CalcFactor(testppxyh, fmd, 1, 1, meas, ARR)
# cfo(res, XYH1_2, veeEuler(wTx1), veeEuler(wTx2))
meas = (vectoarr2(XYH1_2),)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)


# testppxyh(res, fmd, 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end


@testset "test z translation case" begin

##

# dummy value
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

# z translation only
wTx = Vector{AffineMap}(undef,2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(1.0, 0, 0, 0))
wTx[2] = Translation(0, 0, 10.0) ∘ LinearMap(Quat(1.0, 0, 0, 0))

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)


@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

#

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
res = zeros(3)

# cfo = CalcFactor(testppxyh, fmd, 1, 1, meas, ARR)
# cfo(res, XYH1_2, veeEuler(wTx1), veeEuler(wTx2))
meas = (vectoarr2(XYH1_2),)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)


@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10


##

end


@testset "test roll and translate case 1" begin

##

# dummy values
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

# different orientation, roll
wTx = Vector{AffineMap}(undef, 2)
sq2 = 1.0/sqrt(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

# res = zeros(6)
# testpp3(res, fmd, 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

# res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (vectoarr2(XYH1_2),)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)
# res = zeros(3)
# testppxyh(res, fmd, 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))


@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end



@testset "test roll and translate case 2" begin

##

# dummy values
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

wTx = Vector{AffineMap}(undef, 2)
# different orientation, roll
sq2 = 1.0/sqrt(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))
wTx[2] = Translation(0, 10.0, 0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)

# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

# res = zeros(6)
# testpp3(res, fmd, 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))


@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (XYH1_2,)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)
# res = zeros(3)
# testppxyh(res, fmd, 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end


@testset "test pitch and translate case 1" begin

##

# dummy value
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
qq = convert(Rotations.Quat, Rotations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(qq)

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)

# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (XYH1_2,)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)


@show res

@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end



@testset "test pitch and translate case 2" begin

##

# dummy value
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
qq = convert(Rotations.Quat, Rotations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 10.0, 0) ∘ LinearMap(qq)

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (XYH1_2,)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end



@testset "test pitch and translate case 3" begin

##

# dummy values
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)


wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
qq = convert(Rotations.Quat, Rotations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 0, 10.0) ∘ LinearMap(qq)

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

# res = zeros(6)
# testpp3(res, fmd, 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (XYH1_2,)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)
#
# res = zeros(3)
# testppxyh(res, fmd, 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end


@testset "test pitch and translate case 4" begin

##

# dummy values
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
qq = convert(Rotations.Quat, Rotations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 15.0, 10.0) ∘ LinearMap(qq)

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

# res = zeros(6)
# testpp3(res, fmd, 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (XYH1_2,)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)
#
# res = zeros(3)
# testppxyh(res, fmd, 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end


@testset "test yaw and translate case 1" begin

##

# dummy value
fg_ = initfg()
X0 = addVariable!(fg_, :x0, Pose3)
X1 = addVariable!(fg_, :x1, Pose3)

wTx = Vector{AffineMap}(undef, 2)
# different orientation, yaw
qq = convert(Rotations.Quat, Rotations.AngleAxis(pi/2,0,0,1.0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(qq)

# Recalculate XYH
# change toolbox
wTx1 = convert(SE3, wTx[1])
wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level, free orientation
wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
x1Tx2 = (wTx1\wTx2)
wRlx1Tx2 = wRlx1 * x1Tx2
vEx1_2 = veeEuler(wRlx1Tx2)
XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))

meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

# res = zeros(6)
# testpp3(res, fmd, 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = Pose3Pose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

meas = (XYH1_2,)
res = testFactorResidualBinary(testppxyh, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)
#
# res = zeros(3)
# testppxyh(res, fmd, 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))

@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

##

end

end