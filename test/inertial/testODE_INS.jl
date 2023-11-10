# Basic test of DERelative


using Test
using DifferentialEquations
using RoME
import RoME: imuKinematic!

using Interpolations
using IncrementalInference
using Dates
using Statistics
using TensorCast
using StaticArrays
using Manifolds
import Base: convert

##

@testset "DERelative INS Kinematic tests" begin
##

# a user specified ODE in standard form
# inplace `xdot = f(x, u, t)`
# if linear, `xdot = F*x(t) + G*u(t)`
function insKinematic!(dstate, state, u, t)
  # x is a VelPose3 point (assumed ArrayPartition)
  # u is IMU input (assumed [rate; accel])
  Mt = TranslationGroup(3)
  Mr = SpecialOrthogonal(3)
  # convention
  # b is real time body
  # b1 is one discete timestep in the past, equivalent to `r_Sk = r_Sk1 + r_dSk := r_S{k-1} + r_dSk`

  # Using robotics frame fwd-std-dwn <==> North-East-Down
  # ODE cross check taken from Farrell 2008, section 11.2.1, p.388
  # NOTE, Farrell 2008 has gravity logic flipped in some of his equations.
  # WE TAKE gravity as up is positive (think a pendulum hanging back during acceleration)
  # expected gravity in FSD frame (akin to NED).  This is a model of gravity we expect to measure.
  i_G = [0; 0; -9.81] 

  # attitude computer
  w_R_b = state.x[2] # Rotation element
  i_R_b = w_R_b
  # assume body-frame := imu-frame
  b_Ωbi = hat(Mr, Identity(Mr), u[].gyro(t)) # so(3): skew symmetric Lie algebra element
  # assume perfect measurement, i.e. `i` here means measured against native inertial (no coriolis, transport rate, error corrections)
  i_Ṙ_b = i_R_b * b_Ωbi
  # assume world-frame := inertial-frame
  w_Ṙ_b = i_Ṙ_b
  # tie back to the ODE solver

  dstate.x[2] .= w_Ṙ_b
  # Note closed form post integration result (remember exp is solution to first order diff equation)
  # w_R_b = exp(Mr, w_R_b1, b_Ωbi)
  
  # measured inertial acceleration
  b_Abi = u[].accel(t) # already a tangent vector
  # inertial (i.e. world) velocity-dot (accel) by compensating (i.e. removing) expected gravity measurement
  i_V̇ = i_R_b * b_Abi - i_G
  # assume world is inertial frame
  w_V̇ = i_V̇
  dstate.x[3] .= w_V̇ # velocity state
  
  # position-dot (velocity)
  w_V = state.x[3]
  i_V = w_V
  i_Ṗ = i_V
  w_Ṗ = i_Ṗ
  dstate.x[1] .= w_Ṗ
  
  # TODO add biases, see RoME.InertialPose
  # state[4] := gyro bias
  # state[5] := acce bias

  nothing
end

dt = 0.01
N = 101
gyros = [[0.01, 0.0, 0.0] for _ = 1:N]
a0 = [0, 0, -9.81]
accels = [a0]
w_R_b = [1. 0 0; 0 1 0; 0 0 1]
M = SpecialOrthogonal(3)
for g in gyros[1:end-1]
  X = hat(M, Identity(M), g)
  exp!(M, w_R_b, w_R_b, X*dt)
  push!(accels, w_R_b' * a0)
end

gyros_t = linear_interpolation(range(0; step=dt, length=N), gyros)
accels_t = linear_interpolation(range(0; step=dt, length=N), accels)

p = (gyro=gyros_t, accel=accels_t)

p.accel(0.9)

u0 = ArrayPartition([0.0,0,0], Matrix(getPointIdentity(SpecialOrthogonal(3))), [0.,0,0])
tspan = (0.0, 1.0)

prob = ODEProblem(insKinematic!, u0, tspan, Ref(p))

sol = solve(prob)
last(sol)


@test isapprox(last(sol).x[1], [0,0,0]; atol=0.001)
@test isapprox(M, last(sol).x[2], w_R_b; atol=0.001)
@test isapprox(last(sol).x[3], [0,0,0]; atol=0.001)

##
end


@testset "DERelative IMU Kinematic tests" begin
##

# ## TODO consolidate inside module as RoME.imuKinematic
# function imuKinematic!(du, u, p, t)
#   # p is IMU input (assumed [.gyro; .accel])
#   M = SpecialOrthogonal(3)

#   R = u.x[1] # Rotation
#   V = u.x[2] # Velocity 
#   # P = u.x[3] # Position unused here
#   # ω_b = u.x[4] # Gyroscope bias
#   # A_b = u.x[5] # Accelerometer bias

#   ω_m = p[].gyro(t)
#   Ω = hat(M, Identity(M), ω_m)# + ω_b)
#   Ṙ = R * Ω # d/dt R = d/dt exp(Ω*Δt) => Ṙ = exp(ΩΔt)*d(ΩΔt)/dt = exp(ΩΔt)*Ω

#   A_m = p[].accel(t)
#   V̇ = R * (A_m) # + A_b)
#   Ṗ = V

#   du.x[1] .= Ṙ
#   du.x[2] .= V̇
#   du.x[3] .= Ṗ

#   nothing
# end


dt = 0.01
# N = 1001
N = 401
# tspan = (0.0, dt*(N-1))

gyros = [[0.0, 0.0, pi/2] for _ = 1:N]

a0 = [0.0, 0, 0]
accels = [a0]
w_R_b = [1. 0 0; 0 1 0; 0 0 1]
M = SpecialOrthogonal(3)

# b_a = [0.1, 0, 0]
b_a = [0.0, pi/2*10, 0]
for g in gyros[1:end-1]
  X = hat(M, Identity(M), g)
  exp!(M, w_R_b, w_R_b, X*dt)
  push!(accels, b_a .+ w_R_b * a0)
end

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
  imuKinematic!,
  imuForce;
  # state0=allocate(getPointIdentity(RotVelPos)),
  # state1=allocate(getPointIdentity(RotVelPos)),
  dt, 
  problemType=ODEProblem 
);


# cross check on timestamps and tspan used in the ODE problem
@test isapprox.( (9.466848e8, 9.46684804e8), oder.forwardProblem.tspan ) |> all


##

addFactor!( fg, [:w_P0; :w_P1], oder; graphinit=false );



##

@error("WIP testODE_INS.jl")


P1 = approxConvBelief(fg, :w_P0w_P1f1, :w_P1)


# ## basic sample test

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

