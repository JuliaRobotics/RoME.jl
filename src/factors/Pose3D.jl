
# regular prior for Pose2

export PriorPose3, PackedPriorPose3

"""
$(TYPEDEF)

Direct observation information of `Pose3` variable type.
"""
mutable struct PriorPose3{T <: IIF.SamplableBelief} <: IncrementalInference.AbstractPrior
  Zi::T
  # empty constructor
  PriorPose3{T}() where T = new{T}()
  # regular constructor
  PriorPose3{T}(z::T) where {T <: IIF.SamplableBelief} = new{T}(z)
end
# convenience and default object helper
PriorPose3(z::T=MvNormal(zeros(6), LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)]) )) where {T <: SamplableBelief} = PriorPose3{T}(z)

# Standardized sampling function
function getSample(p3::PriorPose3, N::Int=1)
  return (rand(p3.Zi, N),)
end

"""
$(TYPEDEF)

Serialization type for PriorPose3.
"""
mutable struct PackedPriorPose3  <: IncrementalInference.PackedInferenceType
    Zi::String
    PackedPriorPose3() = new()
    PackedPriorPose3(x::AbstractString) = new(x)
end
function convert(::Type{PriorPose3}, packed::PackedPriorPose3)
  # Zi = SE3(d.vecZi[1:3], Quaternion(d.vecZi[4],d.vecZi[5:7]))
  # Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  # return PriorPose3( MvNormal( veeEuler(Zi), Cov) )
  return PriorPose3( extractdistribution(packed.Zi) )
end
function convert(::Type{PackedPriorPose3}, obj::PriorPose3)
  # tf = SE3(d.Zi.μ[1:3], Euler(d.Zi.μ[4:6]...) )
  # v1 = veeQuaternion(tf)
  # v2 = d.Zi.Σ.mat[:];
  # return PackedPriorPose3(v1, v2, size(d.Zi.Σ.mat,1))
  return PackedPriorPose3(string(obj.Zi))
end
