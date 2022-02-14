
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

getManifold(::IIF.InstanceType{<:Pose2Point2BearingRange}) = ProductGroup(ProductManifold(SpecialOrthogonal(2), TranslationGroup(1)))

function getSample(cfo::CalcFactor{<:Pose2Point2BearingRange})
  # defaults, TODO better reuse
  M = getManifold(cfo.factor)
  e0 = ProductRepr([1 0; 0 1.], [0.])

  # vector of tangents
  smpl = hat(M, e0, [rand(cfo.factor.bearing), rand(cfo.factor.range)])

  # return IIF `::Tuple` format
  return smpl
end


function IIF.getMeasurementParametric(s::Pose2Point2BearingRange{<:Normal, <:Normal})

  meas = [mean(s.bearing), mean(s.range)]
  iΣ = [1/var(s.bearing)             0;
                      0  1/var(s.range)]

  return meas, iΣ
end

function (cfo::CalcFactor{<:Pose2Point2BearingRange})(meas, xi, lm)
  SE2 = SpecialEuclidean(2)
  #FIXME fix this factor to work correctly on manifolds
  if false

  Mt = SE2.manifold.manifolds[1]
  Mr = SE2.manifold.manifolds[2]

  # compose rotations
  rRi = xi.parts[2]
  rRo = retract(Mr, rRi, meas.parts[1])

  # new SE2 objects containing the rotation
  rTo = ProductRepr(xi.parts[1], rRo)
  oTl = ProductRepr([meas.parts[2];0], identity_element(Mr, rRi))

  # prediction landmark
  rTo = Manifolds.compose(SE2, rTo, oTl)

  # cartesian difference in predicted and estimated landmark
  δ_l = Manifolds.compose(Mt, inv(Mt, rTo.parts[1]), lm)

  # find the residuals as though in polar coordinates
  δθ = atan(δ_l[2], δ_l[1])
  δr = norm(δ_l)

  return [δθ; δr]

  # FIXME: this is close to the old BR that worked, use to ensure tests are working
  else
    Xi = vee(SE2, xi, log(SE2, identity_element(SE2, xi), xi))

    Mmeas = getManifold(cfo.factor)
    Xmeas = vee(Mmeas, Manifolds.Identity(Mmeas), meas)
    # 1-bearing
    # 2-range
    # world frame
    meas_θ = Xmeas[1]
    meas_r = Xmeas[2]

    θ = meas_θ + Xi[3]
    mx = meas_r*cos(θ)
    my = meas_r*sin(θ)
  
    ex = lm[1] - (mx + Xi[1])
    ey = lm[2] - (my + Xi[2])
  
    # res = [eθ, er]
    eθ = atan((my + Xi[2]), (mx + Xi[1])) - atan(lm[2], lm[1])  # eθ
    er= sqrt(ex^2 + ey^2) # some wasted computation here       # er 
  
    # return hat(facM, ProductRepr([1. 0;0 1], [0.]), [eθ, er])
    return [eθ, er]
  end
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

Base.@kwdef struct PackedPose2Point2BearingRange <: AbstractPackedFactor
    bearstr::PackedSamplableBelief
    rangstr::PackedSamplableBelief
end

function convert( ::Type{PackedPose2Point2BearingRange}, d::Pose2Point2BearingRange )
  return PackedPose2Point2BearingRange( convert(PackedSamplableBelief, d.bearing), convert(PackedSamplableBelief, d.range) )
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert( ::Type{Pose2Point2BearingRange}, d::PackedPose2Point2BearingRange )
  Pose2Point2BearingRange( convert(SamplableBelief, d.bearstr), convert(SamplableBelief, d.rangstr) )
end
