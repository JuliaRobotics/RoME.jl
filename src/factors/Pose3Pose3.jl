# Pose3Pose3 evaluation functions

# ------------------------------------
"""
$(TYPEDEF)

Rigid transform factor between two Pose3 compliant variables.
"""
Base.@kwdef struct Pose3Pose3{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T = MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)]))
end

getManifold(::InstanceType{Pose3Pose3}) = getManifold(Pose3) # Manifolds.SpecialEuclidean(3)

Pose3Pose3(::UniformScaling) = Pose3Pose3()

function (cf::CalcFactor{<:Pose3Pose3})(X, p, q)
    M = getManifold(Pose3)
    q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
    
    # FIXME, should be tangent vector not coordinates -- likely part of ManOpt upgrade
    #TODO allocalte for vee! see Manifolds #412, fix for AD
    # Xc = zeros(6)
    Xc = Vector{eltype(p.x[1])}(undef, 6)
    vee!(M, Xc, q, log(M, q, q̂))
    # Xc = vee(M, q, log(M, q, q̂))
    return Xc
end



#TODO Serialization

"""
$(TYPEDEF)

Serialization type for `Pose3Pose3`.
"""
Base.@kwdef struct PackedPose3Pose3 <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{Pose3Pose3}, packed::PackedPose3Pose3)
  return Pose3Pose3( convert(SamplableBelief, packed.Z) )
end
function convert(::Type{PackedPose3Pose3}, obj::Pose3Pose3)
  return PackedPose3Pose3( convert(PackedSamplableBelief, obj.Z) )
end




#
