
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

function (cfo::CalcFactor{<:Pose2Point2BearingRange})(_measX::AbstractArray{MT}, _p::AbstractArray{PT}, _l::AbstractArray{LT}) where {MT,PT,LT}
  T = promote_type(MT, PT, LT)
  measX = convert(ArrayPartition{T, Tuple{SMatrix{2, 2, T, 4}, SVector{1, T}}}, _measX)
  p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
  l = convert(SVector{2, T}, _l)
  r = cfo(measX, p, l)
  return r
end

function (cfo::CalcFactor{<:Pose2Point2BearingRange})(
                  measX::ArrayPartition{<:Real}, 
                  p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, 
                  l::SVector{2,T}) where T<:Real
  #
  # wl = l
  # wTp = p
  # pl = pTw*wl
  pl  =  transpose(p.x[2]) * (l - p.x[1])
  # δθ = mθ - plθ
  # δr = mr - plr
  δθ = Manifolds.sym_rem(measX.x[1][2] - atan(pl[2], pl[1]))
  δr = measX.x[2][1]  - norm(pl)

  return SA[δθ, δr]
end

# quick check
# pose = (0,0,0),  bear = 0.0,  range = 10.0   ==>  lm = (10,0)
# pose = (0,0,0),  bear = pi/2,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = 0.0,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = pi/2,  range = 10.0   ==>  lm = (-10,0)



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
