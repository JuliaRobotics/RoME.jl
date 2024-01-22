# Basic test of DERelative


using DifferentialEquations
using Interpolations
using Dates
using Statistics
using TensorCast
using StaticArrays
using Manifolds

using IncrementalInference
using RoME
# import RoME: imuKinematic!
# import Base: convert

using Test

##

@testset "DERelative INS Kinematic tests" begin
##

dt = 0.01
N = 101
w_R_b = [1. 0 0; 0 1 0; 0 0 1]
imu = RoME.generateField_InertialMeasurement_RateZ(;
  dt,
  N,
  rate = [0.01, 0, 0],
  accel0 = [0, 0, -9.81], #WHY NEGATIVE, LEGACY, FIX TO +!  DISCREP BETWEEN insKinematic and imuKinematic
  w_R_b
)

gyros_t = linear_interpolation(range(0; step=dt, length=N), imu.gyros)
accels_t = linear_interpolation(range(0; step=dt, length=N), imu.accels)

p = (gyro=gyros_t, accel=accels_t)

p.accel(0.9)

u0 = ArrayPartition([0.0,0,0], Matrix(getPointIdentity(SpecialOrthogonal(3))), [0.,0,0])
tspan = (0.0, 1.0)

prob = ODEProblem(RoME.insKinematic!, u0, tspan, Ref(p))

sol = solve(prob)
last(sol)

M = SpecialOrthogonal(3)
@test isapprox(last(sol).x[1], [0,0,0]; atol=0.001)
@test isapprox(M, last(sol).x[2], w_R_b; atol=0.001)
@test isapprox(last(sol).x[3], [0,0,0]; atol=0.001)

##
end


@testset "DERelative IMU Kinematic tests" begin
##

dt = 0.01
N = 401

imu = RoME.generateField_InertialMeasurement_RateZ(;dt, N)
gyros = imu.gyros
accels = imu.accels

##

fg = initfg()

# the starting points and "0 seconds"
v0 = addVariable!(fg, :w_P0, RotVelPos, timestamp=DateTime(2000,1,1,0,0,0)) 
v1 = addVariable!(fg, :w_P1, RotVelPos, timestamp=DateTime(2000,1,1,0,0,dt*(N-1))) 
# `accurate_time = trunc(getDatetime(var), Second) + (1e-9*getNstime(var) % 1)`

mp = ManifoldPrior(
  getManifold(RotVelPos),
  ArrayPartition(
    SA[ 1.0 0.0 0.0; 
        0.0 1.0 0.0; 
        0.0 0.0 1.0], 
    SA[0.0, 0.0, 0.0], 
    SA[0.0, 0.0, 0.0]
  ),
  MvNormal(diagm(SVector{9}(ones(9)*1e-3)))
)

addFactor!(fg, [:w_P0;], mp)

##


tst = getTimestamp(v0) |> DateTime |> datetime2unix

timestamps = range(tst; step=dt, length=N) # range(0; step=dt, length=N)

gyros_t = linear_interpolation(timestamps, gyros)
accels_t = linear_interpolation(timestamps, accels)

imuForce = Ref((gyro=gyros_t, accel=accels_t))


oder = DERelative(
  fg, 
  [:w_P0; :w_P1], 
  RotVelPos, 
  RoME.imuKinematic!,
  imuForce;
  dt, 
  problemType=ODEProblem 
);


# cross check on timestamps and tspan used in the ODE problem
@test isapprox.( (9.466848e8, 9.46684804e8), oder.forwardProblem.tspan ) |> all


##

addFactor!( fg, [:w_P0; :w_P1], oder; graphinit=false );



##

@error("WIP testODE_INS.jl")

try
  P1 = approxConvBelief(fg, :w_P0w_P1f1, :w_P1)
catch
  @error "FIXME: First try of imuKinematic! convolution failed!"
end

try
  P1 = approxConvBelief(fg, :w_P0w_P1f1, :w_P1)
catch
  @error "FIXME: Second try of imuKinematic! convolution failed!"
end

## basic sample test

# meas = sampleFactor(fg, :x0x1f1, 10)
# @test size(meas[1][1],1) == 1
# @test size(meas,1) == 10


# ## do all forward solutions

# pts = sampleFactor(fg, :x0f1, 100)

# initVariable!(fg, :x0, pts)
# pts_ = approxConv(fg, :x0x1f1, :x1)
# @cast pts[i,j] := pts_[j][i]
# @test 0.3 < Statistics.mean(pts) < 0.4


# ## check that the reverse solve also works

# initVariable!(fg, :x1, pts_)
# pts_ = approxConv(fg, :x0x1f1, :x0)
# @cast pts[i,j] := pts_[j][i]

# # check the reverse solve to be relatively accurate
# ref_ = (getBelief(fg, :x0) |> getPoints)
# @cast ref[i,j] := ref_[j][i]
# @test (pts - ref) |> norm < 1e-4


# ##

# oder_ = DERelative( fg, [:x0; :x3], 
#                     Position{1}, 
#                     firstOrder!,
#                     tstForce, 
#                     dt=0.05, 
#                     problemType=ODEProblem )

# oder_.forwardProblem.u0 .= [1.0]
# sl = DifferentialEquations.solve(oder_.forwardProblem)

# ##


# # Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",layout=(1,1))

# # for lb in [:x0; :x1;:x2;:x3]
# #   x = getTimestamp(getVariable(fg, lb)) |> DateTime |> datetime2unix
# #   xx = [x;x]
# #   yy = [0;1]
# #   Plots.plot!(xx, yy, show=true)
# # end


# ##


# tfg = initfg()
# pts_ = approxConv(fg, :x0f1, :x3, setPPE=true, tfg=tfg)
# # initVariable!(tfg, :x3, pts)


# ##

# @cast pts[i,j] := pts_[j][i]

# @test getPPE(tfg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
# @test getPPE(tfg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
# @test getPPE(tfg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
# @test       Statistics.mean(pts) - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix)[1] < 1.0


# ##

# # plotKDE(tfg, [:x0;:x1;:x2;:x3])


# ## Now test a full solve

# solveTree!(fg);


# ##


# @test getPPE(fg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
# @test getPPE(fg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
# @test getPPE(fg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
# @test getPPE(fg, :x3).suggested - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1


##

end

