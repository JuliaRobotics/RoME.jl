
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
σ_a = 1e-4#0.16e-3*9.81  # noise density m/s²/√Hz
σ_ω = deg2rad(0.0001)  # noise density rad/√Hz
gn = MvNormal(diagm(ones(3)*σ_ω^2 * 1/dt))
an = MvNormal(diagm(ones(3)*σ_a^2 * 1/dt))

# Σy  = diagm([ones(3)*σ_a^2; ones(3)*σ_ω^2])
gyros = [SA[0, 0, 0.001] + rand(gn) for _ = 1:11]
accels = [SA[0, 0, 9.81 - 1] + rand(an) for _ = 1:11]

tst = now(localzone())
tsp = tst + Second(1)
tspan = (tst,tsp)

# a_b = SA[0.,0,0]
# ω_b = SA[0.,0,0]

##

fac = RoME.InertialDynamic(
  tspan,
  dt,
  accels,
  gyros,
)


## build a basic factor graph

fg = initfg()

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


P1 = approxConvBelief(fg, getLabel(f1), :w_P1)


##
end

##