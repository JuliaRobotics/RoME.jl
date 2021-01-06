# test partial pose3 constraints and evaluation

using Statistics
using RoME
using Test

##

N=50
fg = initfg()

v1 = addVariable!(fg,:x1, Pose3, N=N) # 0.001*randn(6,N)
f0 = addFactor!(fg, [:x1;], PriorPose3(MvNormal(zeros(6),1e-2*Matrix{Float64}(LinearAlgebra.I, 6,6))))

sigx = 0.01
sigy = 0.01
sigth = 0.0281
mu1 = [0.0;0.0; -10.0]
prpz = PartialPriorRollPitchZ(
  MvNormal( mu1[1:2], [sigx 0.0; 0.0 sigy]^2 ),
  Normal( mu1[3], sigth )
)

mu2 = [20.0,5.0,pi/2]
xyy = PartialPose3XYYaw(
  MvNormal(
    mu2[1:2],
    [sigx 0.0; 0.0 sigy]^2
  ),
  Normal( mu2[3], sigth )
)

v2 = addVariable!(fg,:x2, Pose3, N=N) # randn(6,N)

f1 = addFactor!(fg, [:x2], prpz, graphinit=false)
f2 = addFactor!(fg, [:x1;:x2], xyy, graphinit=false)

##


@testset "test PartialPriorRollPitchZ evaluations" begin

##

# ensure that at least the first pose is already initialized
doautoinit!(fg, :x1)
@test isInitialized(fg, :x1)

X1pts = getVal(fg, :x1)
@test sum(isnan.(X1pts)) == 0

ppts = approxConv(fg, :x1x2f1,:x2)
@test sum(isnan.(ppts)) == 0


ensureAllInitialized!(fg)
@test isInitialized(fg, :x2)

# get values and ensure that a re-evaluation produces consistent results
X2pts = getVal(fg, :x2)
@test sum(isnan.(X2pts)) == 0

# check that only partial states are updated
pts = IIF.approxConv(fg, :x2f1, :x2, N=N)

newdims = collect(DFG.getSolverData(f1).fnc.usrfnc!.partial)

olddims = setdiff(collect(1:6), newdims)

@test size(pts, 1) == 6
@test size(pts, 2) == N

# check that untouched dimensions (not in the partial list) truely remain untouched
@test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

# check that the prior new dims are updated to new and correct values
# @show Statistics.mean(pts,dims=2)[newdims]
@test sum(abs.(Statistics.mean(pts,dims=2)[newdims]-mu1[[3;1;2]]) .< [0.5; 0.1; 0.1]) == 3

# ensure a forced re-evaluatoin
@test norm(X2pts[newdims,:] - pts[newdims,:]) < 1.0

# memcheck that the exact same values are used
@test norm(X2pts - getVal(v2)) < 1e-10

##

end


@testset "test residual function of PartialPose3XYYaw" begin

##

tfg = initfg()
X0 = addVariable!(tfg, :x0, Pose3)
X1 = addVariable!(tfg, :x1, Pose3)
# fmd = IIF._defaultFactorMetadata([X0;X1])
# res = zeros(3)
# idx = 1
# meas = getSample(xyy)
xi = zeros(6)
xja = zeros(6)
# xyy(res, fmd, idx, meas, xi, xja)

res = testFactorResidualBinary(xyy, Pose3, Pose3, xi, xja)


@test abs(res[1]-mu2[1]) < 0.3
@test abs(res[2]-mu2[2]) < 0.3
@test abs(res[3]-mu2[3]) < 0.2

##

xjb = zeros(6,1)
xjb[collect(xyy.partial),1] = mu2
# res = zeros(3)
# xyy(res, fmd, idx, meas, xi, xjb)

res = testFactorResidualBinary(xyy, Pose3, Pose3, xi, xjb)

@test 0.0 < norm(res) < 0.3

##

addFactor!(tfg, [:x0;:x1], xyy, graphinit=false)

# meas = getSample(xyy,100)
ccw = IIF._getCCW(tfg, :x0x1f1)
meas = freshSamples(ccw, 100)

@test norm(Statistics.std(meas[1],dims=2) - [0.01;0.01;0.002]) < 0.05

##

end



@testset "test PartialPose3XYYaw evaluations" begin

##

# get existing and predict new
X2pts = getBelief(fg, :x2) |> getPoints
pts = approxConv(fg, :x1x2f1, :x2, N=N)

# find which dimensions are and and are not updated by XYYaw partial
newdims = collect(DFG.getSolverData(f2).fnc.usrfnc!.partial)
olddims = setdiff(collect(1:6), newdims)

# check the number of points are correct
@test size(pts, 1) == 6
@test size(pts, 2) == N

# ensure the unchanged dimensions actually remain unchanged
@test norm(X2pts[olddims,:] - pts[olddims,:]) < 1e-10

for i in 1:N
    pts[6,i] = wrapRad(pts[6,i])
end

# ensure the newly updated values match what is specified in mu2
@show mu2 # mu2 is used for XYYaw
# trying to compare world and body frame values -- not right!
@show Statistics.mean(pts[newdims,:],dims=2)
@test sum(abs.(Statistics.mean(pts[newdims,:],dims=2)-mu2) .< [1.5;1.5;0.3]) == 3

# ensure a re-evaluation of the partial factor updates the partial variable dimensions correclty
@test norm(X2pts[newdims,:] - pts[newdims,:]) < 1.0

# ensure that memory pointers are working correctly
memcheck = getVal(v2)
@test norm(X2pts - memcheck) < 1e-10

##

end


@testset "test predictbelief with two functions" begin

##

val, = predictbelief(fg, :x2, ls(fg, :x2), N=N)

for i in 1:N
  val[6,i] = wrapRad(val[6,i])
end

@test size(val, 1) == 6
@test size(val, 2) == N

estmu1mean = Statistics.mean(val[collect(DFG.getSolverData(f1).fnc.usrfnc!.partial),:],dims=2)
estmu2mean = Statistics.mean(val[collect(DFG.getSolverData(f2).fnc.usrfnc!.partial),:],dims=2)

@show estmu1mean
@show estmu2mean
@test sum(abs.(estmu1mean - mu1[[3;1;2]]) .< [0.7; 0.1; 0.1]) == 3
@test sum(abs.(estmu2mean - mu2) .< [0.7; 0.7; 0.15] ) == 3

memcheck = getVal(v2)
@test 1e-10 < norm(val - memcheck)

##

end


#
