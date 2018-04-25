
using RoME
using Base: Test

@testset "Basic Pose2 initialization" begin
  fg = initfg()
  initFactorGraph!(fg, lbl=:x0, ready=0)
  @test true

  N = 100
  fg = initfg()
  initFactorGraph!(fg, init=zeros(3), P0=0.1*eye(3), N=N, lbl=:x0, ready=0,   firstPoseType=Pose2)
  @test true
end


# test Pose3 initialization
@testset "Base Pose3 initialization" begin
  N = 100
  fg = initfg()
  initCov = eye(6)
  [initCov[i,i] = 0.01^2 for i in 4:6];
  initFactorGraph!(fg, P0=initCov, init=zeros(6), N=N, lbl=:x0, ready=0,   firstPoseType=Pose3)
  @test true
end



#
