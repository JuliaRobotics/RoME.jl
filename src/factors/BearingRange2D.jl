
#-------------------------------------------------------------------------------
# bearing and range available

"""
    $TYPEDEF

Bearing and Range constraint from a Pose2 to Point2 variable.
"""
mutable struct Pose2Point2BearingRange{B <: IIF.SamplableBelief, R <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  bearing::B
  range::R
end

getManifold(::Pose2Point2BearingRange) = ProductManifold(SpecialOrthogonal(2), TranslationGroup(1))

function getSample(cfo::CalcFactor{<:Pose2Point2BearingRange}, N::Int=1)

  smpls = [[rand(cfo.factor.bearing), rand(cfo.factor.range)] for _ = 1:N]
  # must return at least first element in `::Tuple` as `::Matrix`
  return (smpls,)
end


function IIF.getMeasurementParametric(s::Pose2Point2BearingRange{<:Normal, <:Normal})

  meas = [mean(s.bearing), mean(s.range)]
  iΣ = [1/var(s.bearing)             0;
                      0  1/var(s.range)]

  return meas, iΣ
end

function (cfo::CalcFactor{<:Pose2Point2BearingRange})(meas, xi, lm)
  SE2 = SpecialEuclidean(2)
  Xi = vee(SE2, xi, log(SE2, identity_element(SE2, xi), xi))

  # 1-bearing
  # 2-range
  # world frame
  θ = meas[1] + Xi[3]
  mx = meas[2]*cos(θ)
  my = meas[2]*sin(θ)

  ex = lm[1] - (mx + Xi[1])
  ey = lm[2] - (my + Xi[2])
  
  # res = [eθ, er]
  eθ = atan((my + Xi[2]), (mx + Xi[1])) - atan(lm[2], lm[1])  # eθ
  er= sqrt(ex^2 + ey^2) # some wasted computation here       # er 

  # return hat(facM, ProductRepr([1. 0;0 1], [0.]), [eθ, er])
  return [eθ, er]
end  
# quick check
# pose = (0,0,0),  bear = 0.0,  range = 10.0   ==>  lm = (10,0)
# pose = (0,0,0),  bear = pi/2,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = 0.0,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = pi/2,  range = 10.0   ==>  lm = (-10,0)


# function (cfo::CalcFactor{<:Pose2Point2BearingRange})(measX, p, q)
#   #
#   M = SpecialEuclidean(2)
#   q_SE = ProductRepr(q, identity_element(SpecialOrthogonal(2), p.parts[2]))

#   X_se2 = log(M, identity_element(M, p), compose(M, inv(M, p), q_SE))
#   X = X_se2.parts[1]
#   # NOTE this: `X̂ = log(M, p, q_SE)` wrong for what we want
#   #go to tangent vector
#   return measX - X 
# end

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
