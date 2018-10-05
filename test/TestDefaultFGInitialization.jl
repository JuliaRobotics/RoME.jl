
using RoME
using Test

@testset "Basic Pose2 initialization" begin
  global fg = initfg()
  initFactorGraph!(fg, lbl=:x0, ready=0)

  global N = 100
  global fg = initfg()
  initFactorGraph!(fg, init=zeros(3), P0=0.1*Matrix{Float64}(LinearAlgebra.I, 3,3), N=N, lbl=:x0, ready=0,   firstPoseType=Pose2)

end


# test Pose3 initialization
@testset "Base Pose3 initialization" begin
  global N = 100
  global fg = initfg()
  global initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
  [initCov[i,i] = 0.01^2 for i in 4:6];
  initFactorGraph!(fg, P0=initCov, init=zeros(6), N=N, lbl=:x0, ready=0, firstPoseType=Pose3)
  @test true
end



#
