
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