
# regular prior for Pose3

"""
$(TYPEDEF)

Direct observation information of `Pose3` variable type.
"""
Base.@kwdef struct PriorPose3{T <: IIF.SamplableBelief, P} <: IncrementalInference.AbstractPrior
  Z::T = MvNormal(zeros(6), diagm([0.01*ones(3);0.0001*ones(3)]))
end

getManifold(::InstanceType{PriorPose3}) = getManifold(Pose3) # SpecialEuclidean(3)


#FIXME Serialization
"""
$(TYPEDEF)

Serialization type for PriorPose3.
"""
Base.@kwdef struct PackedPriorPose3  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{PriorPose3}, packed::PackedPriorPose3)
  return PriorPose3( convert(SamplableBelief, packed.Z) )
end
function convert(::Type{PackedPriorPose3}, obj::PriorPose3)
  return PackedPriorPose3(convert(PackedSamplableBelief, obj.Z))
end
