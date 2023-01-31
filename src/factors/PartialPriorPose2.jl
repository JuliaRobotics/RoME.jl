
"""
$(TYPEDEF)

Constrain only the yaw angle of a Pose2, generally used for Gyrocompass, Magnetometer, Dual-GNSS heading type measurements, or any other similar construct.
"""
Base.@kwdef struct PartialPriorYawPose2{T <: IIF.SamplableBelief} <: IIF.AbstractPrior
  Z::T
  partial::Tuple{Int} = (3,)
end
PartialPriorYawPose2(Z::SamplableBelief) = PartialPriorYawPose2(;Z)

getManifold(::PartialPriorYawPose2) = RealCircleGroup() # SpecialEuclidean(2)

function getSample(cf::CalcFactor{<:PartialPriorYawPose2})
    
  Z = cf.factor.Z
  return rand(Z,1)
  # M = getManifold(cf.factor)
  # p = getPointIdentity(M)
  
  # Xc = [0,0,rand(Z)]
  
  # X = hat(M, p, Xc)
  # points = exp(M, p, X)
  # return points
end

## Serialization support

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPartialPriorYawPose2 <: AbstractPackedFactor
    Z::PackedSamplableBelief
end

function convert(::Type{PackedPartialPriorYawPose2}, d::PartialPriorYawPose2)
  PackedPartialPriorYawPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PartialPriorYawPose2}, d::PackedPartialPriorYawPose2)
  PartialPriorYawPose2(convert(SamplableBelief, d.Z))
end
