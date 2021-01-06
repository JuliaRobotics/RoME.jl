# test PartialPose3XYYaw

using RoME
using CoordinateTransformations, Rotations
# Distributions, TransformUtils
using Test



##

@testset "test x translation case" begin

##

# trivial cases first, orientation based tests below
wTx = Vector{AffineMap}(undef, 2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(1.0, 0, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(Quat(1.0, 0, 0, 0))

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
meas = (veeEuler(x1Tx2),)
res = testFactorResidualBinary(testpp3, Pose3, Pose3, veeEuler(wTx1), veeEuler(wTx2), meas)

@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

#

# test with PartialXYH
testppxyh = PartialPose3XYYaw(
  MvNormal(XYH1_2[1:2], 0.001*diagm([1.0;1.0])),
  Normal(XYH1_2[3], 0.001)
)
# testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*Matrix{Float64}(LinearAlgebra.I, 3,3))  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )

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















#
