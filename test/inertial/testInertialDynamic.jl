
using Test
using RoME
using DifferentialEquations
using Dates, TimeZones
using StaticArrays


##
@testset "test InertialDynamic factor" begin
##


## test factor with rotation around z axis and initial velocity up
# DUPLICATED IN testIMUDeltaFactor.jl

dt = 0.1
N = 11

σ_a = 1e-4 #0.16e-3*9.81  # noise density m/s²/√Hz
σ_ω = deg2rad(0.0001)  # noise density rad/√Hz
imu = RoME.generateField_InertialMeasurement_noise(; dt, N, rate=SA[0, 0, 0.001], accel0=SA[0, 0, 9.81-1], σ_a, σ_ω)

tst = now(localzone())
tsp = tst + Second(imu.tspan[2]-imu.tspan[1])
tspan = (tst,tsp)

##

fac = RoME.InertialDynamic(
  tspan,
  dt,
  imu.accels,
  imu.gyros,
)

## build a basic factor graph

fg = initfg()
getSolverParams(fg).N = 50

addVariable!(fg, :w_P0, RotVelPos; timestamp = tst)
addVariable!(fg, :w_P1, RotVelPos; timestamp = tsp)

addFactor!(fg, [:w_P0;], 
  ManifoldPrior(
    getManifold(RotVelPos),
    ArrayPartition(
      SA[ 1.0 0.0 0.0; 
          0.0 1.0 0.0; 
          0.0 0.0 1.0], 
      SA[10.0, 0.0, 0.0], 
      SA[0.0, 0.0, 0.0]
    ),
    MvNormal(diagm(ones(9)*1e-3))
  )
)

#
f1 = addFactor!(fg, [:w_P0; :w_P1], fac; graphinit=false)

##

@test !isInitialized(fg, :w_P1)
doautoinit!(fg, :w_P0)
@test isInitialized(fg, :w_P0)

##

# flb = getLabel(f1)
# sampleFactor(fg, flb, 50)

##

try
  P1 = approxConvBelief(fg, getLabel(f1), :w_P1)
catch
  @error "FIXME first approxConv on InertialDynamic failed!"
  @test_broken false
end

# P1 = approxConvBelief(fg, getLabel(f1), :w_P1)



##
end

##