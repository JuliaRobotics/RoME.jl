
@kwdef struct InertialDynamic{D<:SamplableBelief} <: AbstractManifoldMinimize
  Z::D
end

getManifold(::InertialDynamic) = getManifold(RotVelPos)





## TODO consolidate inside module as RoME.imuKinematic
function imuKinematic!(du, u, p, t; g=SA[0; 0; 9.81])
  # p is IMU input (assumed [.gyro; .accel])
  M = SpecialOrthogonal(3)

  R = u.x[1] # Rotation
  V = u.x[2] # Velocity 
  # P = u.x[3] # Position unused here
  # ω_b = u.x[4] # Gyroscope bias
  # A_b = u.x[5] # Accelerometer bias

  ω_m = p[].gyro(t)
  Ω = hat(M, Identity(M), ω_m)# + ω_b)
  Ṙ = R * Ω # d/dt R = d/dt exp(Ω*Δt) => Ṙ = exp(ΩΔt)*d(ΩΔt)/dt = exp(ΩΔt)*Ω

  A_m = p[].accel(t)
  V̇ = R * A_m - g # R * (A_m + A_b) - g
  Ṗ = V

  du.x[1] .= Ṙ
  du.x[2] .= V̇
  du.x[3] .= Ṗ

  nothing
end



## ==========================================================
## LEGACY BELOW


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