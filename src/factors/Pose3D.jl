
# regular prior for Pose3

export PriorPose3, PackedPriorPose3

"""
$(TYPEDEF)

Direct observation information of `Pose3` variable type.
"""
struct PriorPose3{T <: IIF.SamplableBelief, P} <: IncrementalInference.AbstractPrior
  Z::T
  p::P
  # empty constructor
  # PriorPose3{T}() where T = new{T}()
  # regular constructor
  # PriorPose3{T}(z::T) where {T <: IIF.SamplableBelief} = new{T}(z)
end
# convenience and default object helper
PriorPose3(z::IIF.SamplableBelief) = PriorPose3(z, getPointIdentity(Pose3))
PriorPose3{T,P}(z::IIF.SamplableBelief) where {T,P} = PriorPose3{T,P}(z, getPointIdentity(Pose3))
PriorPose3() = PriorPose3(MvNormal(zeros(6), LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)])))

DFG.getManifold(::PriorPose3) = SpecialEuclidean(3)

# Standardized sampling function
function getSample(cf::CalcFactor{<:PriorPose3}, N::Int=1)
  Z = cf.factor.Z
  p = cf.factor.p
  M = getManifold(cf.factor)
  
  Xc = [rand(Z) for _ in 1:N]
  
  # X = get_vector.(Ref(M), Ref(p), Xc, Ref(DefaultOrthogonalBasis()))
  X = hat.(Ref(M), Ref(p), Xc)
  points = exp.(Ref(M), Ref(p), X)

  return (points, )
end

#FIXME Serialization
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
  return PriorPose3( convert(SamplableBelief, packed.Zi) )
end
function convert(::Type{PackedPriorPose3}, obj::PriorPose3)
  return PackedPriorPose3(convert(PackedSamplableBelief, obj.Zi))
end
