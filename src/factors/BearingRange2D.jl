
#-------------------------------------------------------------------------------
# bearing and range available

"""
    $TYPEDEF

Bearing and Range constraint from a Pose2 to Point2 variable.
"""
mutable struct Pose2Point2BearingRange{B <: IIF.SamplableBelief, R <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
    bearing::B
    range::R
end

function getSample(cfo::CalcFactor{<:Pose2Point2BearingRange}, N::Int=1)
  smpls = zeros(2, N)
  smpls[1,:] = rand(cfo.factor.bearing, N)[:]
  smpls[2,:] = rand(cfo.factor.range, N)[:]

  # must return at least first element in `::Tuple` as `::Matrix`
  return (smpls,)
end


function IIF.getParametricMeasurement(s::Pose2Point2BearingRange{<:Normal, <:Normal})

  meas = [mean(s.bearing), mean(s.range)]
  iΣ = [1/var(s.bearing)             0;
                      0  1/var(s.range)]

  return meas, iΣ
end

# TODO consolidate with parametric constraint, follow at #467
function (cfo::CalcFactor{<:Pose2Point2BearingRange})(meas, xi, lm)
  # 1-bearing
  # 2-range
  # world frame
  θ = meas[1] + xi[3]
  mx = meas[2]*cos(θ)
  my = meas[2]*sin(θ)

  ex = lm[1] - (mx + xi[1])
  ey = lm[2] - (my + xi[2])
  
  # res = [eθ, er]
  eθ = atan((my + xi[2]), (mx + xi[1])) - atan(lm[2], lm[1])  # eθ
  er= sqrt(ex^2 + ey^2) # some wasted computation here       # er 

  # IIF v0.21+
  return [eθ, er]
end
# quick check
# pose = (0,0,0),  bear = 0.0,  range = 10.0   ==>  lm = (10,0)
# pose = (0,0,0),  bear = pi/2,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = 0.0,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = pi/2,  range = 10.0   ==>  lm = (-10,0)

# Support for database based solving

passTypeThrough(d::FunctionNodeData{Pose2Point2Range}) = d

mutable struct PackedPose2Point2BearingRange <: IncrementalInference.PackedInferenceType
    bearstr::String
    rangstr::String
end

function convert( ::Type{PackedPose2Point2BearingRange}, d::Pose2Point2BearingRange )
  return PackedPose2Point2BearingRange( convert(PackedSamplableBelief, d.bearing), convert(PackedSamplableBelief, d.range) )
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert( ::Type{Pose2Point2BearingRange}, d::PackedPose2Point2BearingRange )
  Pose2Point2BearingRange( convert(SamplableBelief, d.bearstr), convert(SamplableBelief, d.rangstr) )
end
