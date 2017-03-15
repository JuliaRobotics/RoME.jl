# test partial pose3 constraints and evaluation

using RoME, Distributions
using IncrementalInference, TransformUtils
using Base: Test

N=50
fg = initfg()

v1 = addNode!(fg,:x1,0.001*randn(6,N),N=N)

mu1 = [2.0;0.0; -10.0]
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

v2 = addNode!(fg,:x2,randn(6,N),N=N)

f1 = addFactor!(fg, [:x2], prpz)
f2 = addFactor!(fg, [:x1;:x2], xyy)



println("test PartialPriorRollPitchZ evaluations")

X2pts = getVal(v2)
pts = evalFactor2(fg, f1, v2.index, N=N)

newdims = collect(getData(f1).fnc.usrfnc!.partial)
olddims = setdiff(collect(1:6), newdims)

@test size(pts, 1) == 6
@test size(pts, 2) == N
@test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10
# @show Base.mean(pts,2)[newdims]
@test norm(Base.mean(pts,2)[newdims]-mu1) < 0.3
# ensure the correct response from
@test norm(X2pts[newdims,:] - pts[newdims,:]) > 2.0
# memcheck
@test norm(X2pts - getVal(v2)) < 1e-10



println("test residual function of PartialPose3XYYaw")

res = zeros(3)
idx = 1
meas = getSample(xyy)
xi = zeros(6,1)
xja = zeros(6,1)
xyy(res, idx, meas, xi, xja)
@test abs(res[1]-mu2[1]) < 0.2
@test abs(res[2]-mu2[2]) < 0.2
@test abs(res[3]-mu2[3]) < 0.2

xjb = zeros(6,1)
xjb[collect(xyy.partial),1] = mu2
res = zeros(3)
xyy(res, idx, meas, xi, xjb)
@test 0.0 < norm(res) < 0.2

meas = getSample(xyy,100)
@test norm(Base.std(meas[1],2)- [0.01;0.01;0.002]) < 0.005




println("test PartialPose3XYYaw evaluations")

X2pts = getVal(v2)
pts = evalFactor2(fg, f2, v2.index, N=N)

newdims = collect(getData(f2).fnc.usrfnc!.partial)
olddims = setdiff(collect(1:6), newdims)

@test size(pts, 1) == 6
@test size(pts, 2) == N
@test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

# @show Base.mean(getVal(v1),2)[newdims]
# @show mu2
# @show Base.mean(pts,2)[newdims]
# @show Base.std(pts,2)[newdims]
for i in 1:N
  pts[6,i] = wrapRad(pts[6,i])
end

@test norm(Base.mean(pts[newdims,:],2)-mu2) < 0.3
# ensure the correct response from
@test norm(X2pts[newdims,:] - pts[newdims,:]) > 2.0
memcheck = getVal(v2)
@test norm(X2pts - memcheck) < 1e-10




println("test predictbelief with two functions")
val = predictbelief(fg, :x2, [:x2;:x1x2], N=N)

for i in 1:N
  val[6,i] = wrapRad(val[6,i])
end

@test size(val, 1) == 6
@test size(val, 2) == N


@test norm(Base.mean(val[collect(getData(f1).fnc.usrfnc!.partial),:],2)-mu1) < 0.3
@test norm(Base.mean(val[collect(getData(f2).fnc.usrfnc!.partial),:],2)-mu2) < 0.3

memcheck = getVal(v2)
@test 1e-10 < norm(val - memcheck)



#
