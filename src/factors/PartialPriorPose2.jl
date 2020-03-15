
"""
$(TYPEDEF)

Constrain only the yaw angle of a Pose2, generally used for Gyrocompass, Magnetometer, Dual-GNSS heading type measurements, or any other similar construct.
"""
mutable struct PartialPriorYawPose2{T} <: IIF.FunctorSingleton  where {T <: IIF.SamplableBelief}
    Z::T
    partial::Tuple{Int}
    PartialPriorYawPose2{T}() where T = new{T}()
    PartialPriorYawPose2{T}(x::T) where {T <: IIF.SamplableBelief}  = new{T}(x, (3,))
end
PartialPriorYawPose2(x::T) where {T <: IIF.SamplableBelief} = PartialPriorYawPose2{T}(x)

function getSample(p2::PartialPriorYawPose2, N::Int=1)
  return (reshape(rand(p2.Z,N),1,N), )
end



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
