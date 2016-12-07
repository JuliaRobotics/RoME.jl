# test camera functions

using RoME, IncrementalInference
using Base.Test

ci = CameraIntrinsic()
ce = CameraExtrinsic()
pt = [1.0;0.0;5.0]

tol = 1e-8

gg = (x, res) -> cameraResidual!(res, x, ci, ce, pt)
# res = zeros(2)
# @time gg([0.0;0.0], res)

# Profile.clear()
# @profile
y = numericRootGenericRandomizedFnc(
        gg,
        2, 2, randn(2)  )
#
@test abs((ci.K[1,3]+ci.K[1,1]*pt[1]/pt[3]) - y[1]) < tol
@test abs(y[2] - ci.K[2,3]) < tol



# using cameraResidual function
gg = (x, res) -> cameraResidual!(res, x, ci, ce, pt)
# res = zeros(2)
# @time gg([0.0;0.0], res)

# Profile.clear()
# @profile
y = numericRootGenericRandomizedFnc(
        gg,
        2, 2, randn(2)  )
#
@test abs((ci.K[1,3]+ci.K[1,1]*pt[1]/pt[3]) - y[1]) < tol
@test abs(y[2] - ci.K[2,3]) < tol

















#
