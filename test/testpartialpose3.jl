# test partial pose3 constraints and evaluation

using RoME, Distributions
using IncrementalInference
using Base: Test

N=50
fg = initfg()

v1 = addNode!(fg,:x1,0.01*randn(6,N),N=N)

mu1 = [0.0173329;0.158899; -10.7]
prpz = PartialPriorRollPitchZ(
  MvNormal( mu1[1:2], [[2.8092e-6;0.0]'; [0.0;2.8092e-6]'] ),
  Normal(mu1[3], 0.0281)
)

mu2 = [-1.72616,-0.126154,0.0380778]
xyy = PartialPose3XYYaw(
  MvNormal(
    mu2,
    [0.150287 0.0 0.0; 0.0 0.150287 0.0; 0.0 0.0 6.39218e-8]^2
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

@test norm(Base.mean(pts,2)[newdims]-mu2) < 0.3
# ensure the correct response from
@test norm(X2pts[newdims,:] - pts[newdims,:]) > 2.0
memcheck = getVal(v2)
@test norm(X2pts - memcheck) < 1e-10















#
