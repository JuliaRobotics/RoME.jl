module RoMEDifferentialEquationsExt

using DifferentialEquations
using Manifolds
using Dates, TimeZones
using Interpolations
using TensorCast
# using RoME
import RoME: imuKinematic!, RotVelPos, InertialDynamic
import RoME: getPointIdentity
import IncrementalInference: DERelative

function InertialDynamic(
  tspan::Tuple{<:Real, <:Real}, # = _calcTimespan(Xi),
  dt::Real,
  gyros::AbstractVector,
  accels::AbstractVector;
  N::Integer = size(gyros,1),
  timestamps = collect(range(tspan[1]; step=dt, length=N)),
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

# function InertialDynamic(
#   tspan::Tuple{<:Real, <:Real}, # = _calcTimespan(Xi),
#   dt::Real,
#   gyros::AbstractVector,
#   accels::AbstractVector;
#   kw...
# )
#   @cast gyros_[i,d] := gyros[i][d]
#   @cast accels_[i,d] := accels[i][d]
#   InertialDynamic(tspan,dt,gyros_,accels_;kw...)
# end

function InertialDynamic(
  tspan::Tuple{<:ZonedDateTime, <:ZonedDateTime},
  w...;
  kw...
)
  tspan_ = map(t -> datetime2unix(DateTime(t)), tspan)
  InertialDynamic(tspan_, w...; kw...)
end


end # weakmod