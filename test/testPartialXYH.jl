# test PartialPose3XYYaw

using RoME, CoordinateTransformations, Distributions, TransformUtils
using Base: Test

println("test x translation case")
# trivial cases first, orientation based tests below
wTx = Vector{AffineMap}(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(1.0, 0, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(Quat(1.0, 0, 0, 0))

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10







println("test z translation case")
# z translation only
wTx = Vector{AffineMap}(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(1.0, 0, 0, 0))
wTx[2] = Translation(0, 0, 10.0) ∘ LinearMap(Quat(1.0, 0, 0, 0))

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10







println("test roll and translate case 1")
# different orientation, roll
sq2 = 1.0/sqrt(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10








println("test roll and translate case 2")
# different orientation, roll
sq2 = 1.0/sqrt(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))
wTx[2] = Translation(0, 10.0, 0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10










println("test pitch and translate case 1")
# different orientation, pitch
qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(qq)

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10








println("test pitch and translate case 2")
# different orientation, pitch
qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 10.0, 0) ∘ LinearMap(qq)

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10










println("test pitch and translate case 3")
# different orientation, pitch
qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 0, 10.0) ∘ LinearMap(qq)

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10











println("test pitch and translate case 4")
# different orientation, pitch
qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 15.0, 10.0) ∘ LinearMap(qq)

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10















println("test yaw and translate case 1")
# different orientation, yaw
qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/2,0,0,1.0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(qq)

## Recalculate XYH
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
testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*eye(6)))
res = zeros(6)
testpp3(res, 1, (veeEuler(x1Tx2),), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*eye(3))  )
res = zeros(3)
testppxyh(res, 1, ((XYH1_2')',), (veeEuler(wTx1)')', (veeEuler(wTx2)')')
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

















#
