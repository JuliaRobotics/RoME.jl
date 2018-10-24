# test PartialPose3XYYaw

using RoME
using CoordinateTransformations
# Distributions, TransformUtils
using Test

@testset "test x translation case" begin

# trivial cases first, orientation based tests below
global wTx = Vector{AffineMap}(undef, 2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(1.0, 0, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(Quat(1.0, 0, 0, 0))

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(
  MvNormal(XYH1_2[1:2], 0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)),
  Normal(XYH1_2[3], 0.001)
)
# global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2,0.001*Matrix{Float64}(LinearAlgebra.I, 3,3))  )

global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end





@testset "test z translation case" begin

# z translation only
global wTx = Vector{AffineMap}(undef,2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(1.0, 0, 0, 0))
wTx[2] = Translation(0, 0, 10.0) ∘ LinearMap(Quat(1.0, 0, 0, 0))

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end





@testset "test roll and translate case 1" begin

# different orientation, roll
global wTx = Vector{AffineMap}(undef, 2)
global sq2 = 1.0/sqrt(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end






@testset "test roll and translate case 2" begin

global wTx = Vector{AffineMap}(undef, 2)
# different orientation, roll
global sq2 = 1.0/sqrt(2)
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))
wTx[2] = Translation(0, 10.0, 0) ∘ LinearMap(Quat(sq2, sq2, 0, 0))

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end








@testset "test pitch and translate case 1" begin

global wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
global qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(qq)

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end






@testset "test pitch and translate case 2" begin

global wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
global qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 10.0, 0) ∘ LinearMap(qq)

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end








@testset "test pitch and translate case 3" begin

global wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
global qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 0, 10.0) ∘ LinearMap(qq)

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end









@testset "test pitch and translate case 4" begin

global wTx = Vector{AffineMap}(undef, 2)
# different orientation, pitch
global qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/4,0,1.0,0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(0, 15.0, 10.0) ∘ LinearMap(qq)

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end













@testset "test yaw and translate case 1" begin

global wTx = Vector{AffineMap}(undef, 2)
# different orientation, yaw
global qq = convert(CoordinateTransformations.Quat, CoordinateTransformations.AngleAxis(pi/2,0,0,1.0))
wTx[1] = Translation(0.0,0,0) ∘ LinearMap(qq)
wTx[2] = Translation(10.0, 0, 0) ∘ LinearMap(qq)

## Recalculate XYH
# change toolbox
global wTx1 = convert(SE3, wTx[1])
global wTx2 = convert(SE3, wTx[2])
# get rotation to world frame for local level, free orientation
global wEx1 = convert(Euler, wTx1.R)
wEx1.Y = 0.0
global wRlx1 = SE3(zeros(3), wEx1)
# wRx2 = deepcopy(wTx2); wRx2.t = zeros(3);
# Odometries
global x1Tx2 = (wTx1\wTx2)
global wRlx1Tx2 = wRlx1 * x1Tx2
global vEx1_2 = veeEuler(wRlx1Tx2)
global XYH1_2 = [vEx1_2[1], vEx1_2[2], vEx1_2[6]]
global prz2 = veeEuler(wTx1)[[3;4;5]]

# test with Pose3Pose3
global testpp3 = Pose3Pose3(MvNormal(veeEuler(x1Tx2),0.001*Matrix{Float64}(LinearAlgebra.I, 6,6)))
global res = zeros(6)
testpp3(res, FactorMetadata(), 1, (veeEuler(x1Tx2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:3]) < 1e-10
@test norm(res[4:6]) < 1e-10

# test with PartialXYH
global testppxyh = PartialPose3XYYaw(  MvNormal(XYH1_2[1:2],0.001*Matrix{Float64}(LinearAlgebra.I, 2,2)), Normal(XYH1_2[3], 0.001)  )
global res = zeros(3)
testppxyh(res, FactorMetadata(), 1, (vectoarr2(XYH1_2),), vectoarr2(veeEuler(wTx1)), vectoarr2(veeEuler(wTx2)))
@show res
@test norm(res[1:2]) < 1e-10
@test abs(res[3]) < 1e-10

end















#
