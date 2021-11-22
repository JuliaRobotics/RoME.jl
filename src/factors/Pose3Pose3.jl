# Pose3Pose3 evaluation functions

# ------------------------------------
"""
$(TYPEDEF)

Rigid transform factor between two Pose3 compliant variables.
"""
struct Pose3Pose3{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T
end
# convenience and default constructor
Pose3Pose3() = Pose3Pose3(MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)])))

DFG.getManifold(::Pose3Pose3) = Manifolds.SpecialEuclidean(3)

Pose3Pose3(::UniformScaling) = Pose3Pose3()

function (cf::CalcFactor{<:Pose3Pose3})(X, p, q)
    M = getManifold(Pose3)
    q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
    #TODO allocalte for vee! see Manifolds #412, fix for AD
    Xc = zeros(6)
    vee!(M, Xc, q, log(M, q, q̂))
    return Xc
end



#TODO Serialization

"""
$(TYPEDEF)

Serialization type for `Pose3Pose3`.
"""
mutable struct PackedPose3Pose3 <: IncrementalInference.PackedInferenceType
  Z::String
  PackedPose3Pose3() = new()
  PackedPose3Pose3(x::AbstractString) = new(x)
end
function convert(::Type{Pose3Pose3}, packed::PackedPose3Pose3)
  return Pose3Pose3( convert(SamplableBelief, packed.Z) )
end
function convert(::Type{PackedPose3Pose3}, obj::Pose3Pose3)
  return PackedPose3Pose3( convert(PackedSamplableBelief, obj.Z) )
end




#
