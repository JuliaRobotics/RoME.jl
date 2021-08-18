
"""
$(TYPEDEF)

Constrain only the yaw angle of a Pose2, generally used for Gyrocompass, Magnetometer, Dual-GNSS heading type measurements, or any other similar construct.
"""
mutable struct PartialPriorYawPose2{T <: IIF.SamplableBelief} <: IIF.AbstractPrior
    Z::T
    partial::Tuple{Int}
    PartialPriorYawPose2{T}() where T = new{T}()
    PartialPriorYawPose2{T}(x::T) where {T <: IIF.SamplableBelief}  = new{T}(x, (3,))
end
PartialPriorYawPose2(x::T) where {T <: IIF.SamplableBelief} = PartialPriorYawPose2{T}(x)

function getSample(cf::CalcFactor{<:PartialPriorYawPose2}, N::Int=1)
    
  Z = cf.factor.Z
  M = getManifold(cf.factor)
  p = identity_element(M)
  
  Xc = [rand(Z) for _ in 1:N]
  
  X = hat.(Ref(M), Ref(p), Xc)
  points = exp.(Ref(M), Ref(p), X)
  return (points, )
end

getManifold(::PartialPriorYawPose2) = SpecialOrthogonal(2)

## Serialization support

"""
$(TYPEDEF)
"""
mutable struct PackedPartialPriorYawPose2 <: IIF.PackedInferenceType
    Z::String
    PackedPartialPriorYawPose2() = new()
    PackedPartialPriorYawPose2(x::T) where {T <: AbstractString}  = new(x)
end

function convert(::Type{PackedPartialPriorYawPose2}, d::PartialPriorYawPose2)
  PackedPartialPriorYawPose2(string(d.Z))
end
function convert(::Type{PartialPriorYawPose2}, d::PackedPartialPriorYawPose2)
  PartialPriorYawPose2(extractdistribution(d.Z))
end
