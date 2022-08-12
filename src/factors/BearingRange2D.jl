
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
  e0 = ArrayPartition([1 0; 0 1.], [0.])

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

function (cfo::CalcFactor{<:Pose2Point2BearingRange})(measX, p, l)
  #
  M = getManifold(cfo.factor)

  # wl = l
  # wTp = p
  # pl = pTw*wl
  pl  =  transpose(p.x[2]) * (l - p.x[1])

  mθ,mr = vee(M, Manifolds.Identity(M), measX) 
  δθ = mθ - atan(pl[2], pl[1])
  δr = mr - norm(pl)

  return [δθ; δr]
end

# quick check
# pose = (0,0,0),  bear = 0.0,  range = 10.0   ==>  lm = (10,0)
# pose = (0,0,0),  bear = pi/2,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = 0.0,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = pi/2,  range = 10.0   ==>  lm = (-10,0)



# Support for database based solving

passTypeThrough(d::FunctionNodeData{<:Pose2Point2Range}) = d

Base.@kwdef struct PackedPose2Point2BearingRange <: AbstractPackedFactor
    bearstr::PackedSamplableBelief
    rangstr::PackedSamplableBelief
end

function convert( ::Type{<:PackedPose2Point2BearingRange}, d::Pose2Point2BearingRange )
  return PackedPose2Point2BearingRange( convert(PackedSamplableBelief, d.bearing), convert(PackedSamplableBelief, d.range) )
end

function convert( ::Type{<:Pose2Point2BearingRange}, d::PackedPose2Point2BearingRange )
  Pose2Point2BearingRange( convert(SamplableBelief, d.bearstr), convert(SamplableBelief, d.rangstr) )
end
