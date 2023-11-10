module RoMEDifferentialEquationsExt

using DifferentialEquations
using Manifolds
using Interpolations
# using RoME
import RoME: imuKinematic!, RotVelPos, InertialDynamic
import RoME: getPointIdentity

function InertialDynamic(
  tspan::Tuple{<:Real, <:Real}, # = _calcTimespan(Xi),
  dt::Real,
  gyros::AbstractMatrix,
  accels::AbstractMatrix;
  N::Integer = size(gyros,1),
  timestamps = range(tspan[1]; step=dt, length=N) # range(0; step=dt, length=N)
)
  # use data interpolation
  gyros_t = linear_interpolation(timestamps, gyros)
  accels_t = linear_interpolation(timestamps, accels)  

  data = Ref((gyro=gyros_t, accel=accels_t))

  domain = RotVelPos
  state0 = allocate(getPointIdentity(domain))
  state1 = allocate(getPointIdentity(domain))
  
  problemType = ODEProblem
  # forward time problem
  fproblem = problemType(imuKinematic!, state0, tspan, data; dt)
  # backward time problem
  bproblem = problemType(imuKinematic!, state1, (tspan[2], tspan[1]), data; dt = -dt)

  # build the IIF recognizable object
  return DERelative(domain, fproblem, bproblem, data)
end


end # weakmod