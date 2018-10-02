# test partial pose3 constraints and evaluation

using RoME, Distributions
using IncrementalInference, TransformUtils
using Base: Test

N=50
fg = initfg()

v1 = addNode!(fg,:x1, Pose3, N=N) # 0.001*randn(6,N)
f0 = addFactor!(fg, [:x1;], PriorPose3(MvNormal(zeros(6),1e-2*eye(6))))

mu1 = [0.0;0.0; -10.0]
prpz = PartialPriorRollPitchZ(
  MvNormal( mu1[1:2], [[1e-4;0.0]'; [0.0;1e-4]'] ),
  Normal( mu1[3], 0.0281 )
)

mu2 = [20.0,5.0,pi/2]
xyy = PartialPose3XYYaw(
  MvNormal(
    mu2,
    [1e-4 0.0 0.0; 0.0 1e-4 0.0; 0.0 0.0 4e-6]
  )
)

v2 = addNode!(fg,:x2, Pose3, N=N) # randn(6,N)

f1 = addFactor!(fg, [:x2], prpz)
f2 = addFactor!(fg, [:x1;:x2], xyy)


# ls(fg, :x2)

@testset "test PartialPriorRollPitchZ evaluations" begin

  # ensure that at least the first pose is already initialized
  @test isInitialized(fg, :x1)
  ensureAllInitialized!(fg)
  @test isInitialized(fg, :x2)

  # get values and ensure that a re-evaluation produces consistent results
  X2pts = getVal(v2)
  pts = evalFactor2(fg, f1, v2.index, N=N)

  newdims = collect(getData(f1).fnc.usrfnc!.partial)
  olddims = setdiff(collect(1:6), newdims)

  @test size(pts, 1) == 6
  @test size(pts, 2) == N

  # check that untouched dimensions (not in the partial list) truely remain untouched
  @test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

  # check that the prior new dims are updated to new and correct values
  # @show Base.mean(pts,2)[newdims]
  @test sum(abs.(Base.mean(pts,2)[newdims]-mu1[[3;1;2]]) .< [0.5; 0.1; 0.1]) == 3

  # ensure a forced re-evaluatoin
  @test norm(X2pts[newdims,:] - pts[newdims,:]) < 1.0

  # memcheck that the exact same values are used
  @test norm(X2pts - getVal(v2)) < 1e-10
end


@testset "test residual function of PartialPose3XYYaw" begin

  res = zeros(3)
  idx = 1
  meas = getSample(xyy)
  xi = zeros(6,1)
  xja = zeros(6,1)
  xyy(res, nothing, idx, meas, xi, xja)
  @test abs(res[1]-mu2[1]) < 0.2
  @test abs(res[2]-mu2[2]) < 0.2
  @test abs(res[3]-mu2[3]) < 0.2

  xjb = zeros(6,1)
  xjb[collect(xyy.partial),1] = mu2
  res = zeros(3)
  xyy(res, nothing, idx, meas, xi, xjb)
  @test 0.0 < norm(res) < 0.2

  meas = getSample(xyy,100)
  @test norm(Base.std(meas[1],2) - [0.01;0.01;0.002]) < 0.005
end



@testset "test PartialPose3XYYaw evaluations" begin

  X2pts = getVal(v2)
  pts = evalFactor2(fg, f2, v2.index, N=N)

  newdims = collect(getData(f2).fnc.usrfnc!.partial)
  olddims = setdiff(collect(1:6), newdims)

  @test size(pts, 1) == 6
  @test size(pts, 2) == N

  # ensure the unchanged dimensions actually remain unchanged
  @test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

  # @show Base.mean(getVal(v1),2)[newdims]
  # @show mu2
  # @show Base.mean(pts,2)[newdims]
  # @show Base.std(pts,2)[newdims]
  for i in 1:N
    pts[6,i] = wrapRad(pts[6,i])
  end

  # ensure the newly updated values match what is specified in mu2
  @show abs.(Base.mean(pts[newdims,:],2)-mu2)
  @test sum(abs.(Base.mean(pts[newdims,:],2)-mu2) .< [0.7;0.7;0.2]) == 3

  # ensure a re-evaluation of the partial factor updates the partial variable dimensions correclty
  @test norm(X2pts[newdims,:] - pts[newdims,:]) < 1.0

  # ensure that memory pointers are working correctly
  memcheck = getVal(v2)
  @test norm(X2pts - memcheck) < 1e-10
end


@testset "test predictbelief with two functions" begin
  val = predictbelief(fg, :x2, ls(fg, :x2), N=N)

  for i in 1:N
    val[6,i] = wrapRad(val[6,i])
  end

  @test size(val, 1) == 6
  @test size(val, 2) == N

  estmu1mean = Base.mean(val[collect(getData(f1).fnc.usrfnc!.partial),:],2)
  estmu2mean = Base.mean(val[collect(getData(f2).fnc.usrfnc!.partial),:],2)

  @test sum(abs.(estmu1mean - mu1[[3;1;2]]) .< [0.7; 0.1; 0.1]) == 3
  @test sum(abs.(estmu2mean - mu2) .< [0.7; 0.7; 0.15] ) == 3

  memcheck = getVal(v2)
  @test 1e-10 < norm(val - memcheck)
end


#
