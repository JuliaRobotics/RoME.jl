

# using Revise
using Test
using CameraModels
using RoME
using Manifolds

# using ManifoldDiff
# import FiniteDifferences as FD


##
@testset "GenericProjection of 2 camera poses to a point" begin
##

cam = CameraModels.CameraCalibration()

obs1 = CameraModels.PixelIndex(240, 320)
obs2 = CameraModels.PixelIndex(240, 315)

w_T_c1 = ArrayPartition([0; 0  ;0.],[0 0 1; -1 0 0; 0 -1 0.])
w_T_c2 = ArrayPartition([0;-0.1;0.],[0 0 1; -1 0 0; 0 -1 0.])

# w
# [
# 0 -0.1 0
# ]
# =
# w
# [
#  0  0  1
# -1  0  0
#  0 -1  0
# ]
# c
# [
# 0.1 0 0
# ]

##

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :w_P_c1, Pose3)
addVariable!(fg, :w_P_c2, Pose3)
addVariable!(fg, :w_Ph, Point3)

c1 = manikde!(Pose3, [w_T_c1 for _ in 1:1], bw=ones(6));
c2 = manikde!(Pose3, [w_T_c2 for _ in 1:1], bw=ones(6));
h1 = manikde!(Position3, [zeros(3) for _ in 1:1], bw=ones(3));

initVariable!(fg, :w_P_c1, c1, :parametric)
initVariable!(fg, :w_P_c2, c2, :parametric)
initVariable!(fg, :w_Ph, h1, :parametric)

Z=MvNormal([240.0;320], [1 0; 0 1.])
f1 = RoME.GenericProjection{Pose3,Point3}(cam, Z)
addFactor!(fg, [:w_P_c1; :w_Ph], f1)

Z=MvNormal([240.0;315], [1 0; 0 1.])
f2 = RoME.GenericProjection{Pose3,Point3}(cam, Z)
addFactor!(fg, [:w_P_c2; :w_Ph], f2)

# FIXME, should be done on parametric initVariable above
getSolverData(fg[:w_P_c1], :parametric).bw = diagm(ones(6))
getSolverData(fg[:w_P_c2], :parametric).bw = diagm(ones(6))
getSolverData(fg[:w_Ph], :parametric).bw = diagm(ones(3))

M = getManifold(fg, :w_P_c1)
addFactor!(
  fg, 
  [:w_P_c1,], 
  PriorPose3(
    MvNormal(
      vee(M, identity_element(M), log(M,identity_element(M),w_T_c1)),
      diagm(0.01*ones(6))
    )
  )
)

addFactor!(
  fg, 
  [:w_P_c2,], 
  PriorPose3(
    MvNormal(
      vee(M, identity_element(M), log(M,identity_element(M),w_T_c2)),
      diagm(0.01*ones(6))
    )
  )
)

@error "TODO Work in progress for solving GenericProjection factors via solveGraphParametric!"
# IIF.solveGraphParametric!(fg)

##

# using test development function
w_P3 = solveMultiviewLandmark!(fg, :w_Ph)


##


@test isapprox([10.56;0;0], w_P3; atol=1e-3)


##

filepath = joinpath(tempdir(), "testgeneric.tar.gz")
saveDFG(filepath, fg)
fg_ = loadDFG!(initfg(), filepath)

Base.rm(filepath)

##

w_P3 = solveMultiviewLandmark!(fg_, :w_Ph)
@test isapprox([10.56;0;0], w_P3; atol=1e-3)

##
end


##