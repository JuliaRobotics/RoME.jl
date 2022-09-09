
struct P2P2BearingReuse
  measvec::Vector{Float64}
  predvec::Vector{Float64}
  resid::Vector{Float64}
  P2P2BearingReuse() = new(zeros(2),zeros(2),zeros(2))
end

"""
    $TYPEDEF

Single dimension bearing constraint from Pose2 to Point2 variable.
"""
Base.@kwdef struct Pose2Point2Bearing{B <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
    Z::B = Normal()
end

preambleCache(::AbstractDFG, ::AbstractVector{<:DFGVariable}, ::Pose2Point2Bearing) = P2P2BearingReuse()

getManifold(::Pose2Point2Bearing) = SpecialOrthogonal(2)
# Pose2Point2Bearing(x1::B) where {B <: IIF.SamplableBelief} = Pose2Point2Bearing{B}(x1)
# FIXME, there might be something wrong with the sign here
function getSample(cfo::CalcFactor{<:Pose2Point2Bearing})
  return rand(cfo.factor.Z)
end

function (cfo::CalcFactor{<:Pose2Point2Bearing})(Xc, p, l)
  # wl = l
  # wTp = p
  # pl = pTw*wl
  pl  =  transpose(p.x[2]) * (l - p.x[1])
  # δθ = mθ - plθ
  δθ = Manifolds.sym_rem(Xc[1] - atan(pl[2], pl[1]))
  return [δθ]
end


# Packing and Unpacking
Base.@kwdef struct PackedPose2Point2Bearing <: AbstractPackedFactor
    Z::PackedSamplableBelief
    # _type::String = "RoME.PackedPose2Point2Bearing"
end
function convert(::Type{PackedPose2Point2Bearing}, d::Pose2Point2Bearing{B}) where {B <: IIF.SamplableBelief}
  return PackedPose2Point2Bearing(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Pose2Point2Bearing}, d::PackedPose2Point2Bearing)
  Pose2Point2Bearing(convert(SamplableBelief, d.Z))
end
