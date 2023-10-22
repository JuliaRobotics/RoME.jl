


"""
$(TYPEDEF)
"""
mutable struct Point2Point2Range{D <: IIF.SamplableBelief} <: IncrementalInference.AbstractManifoldMinimize # AbstractRelativeMinimize
  Z::D
end

getManifold(::InstanceType{Point2Point2Range}) = TranslationGroup(1)


function (cfo::CalcFactor{<:Point2Point2Range})(rho, xi, lm)
  # Basically `EuclidDistance`
  # must return all dimensions
  return rho .- norm(lm[1:2] .- xi[1:2])
end

passTypeThrough(d::FunctionNodeData{Point2Point2Range}) = d

"""
$(TYPEDEF)
"""
Base.@kwdef mutable struct PackedPoint2Point2Range  <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{PackedPoint2Point2Range}, d::Point2Point2Range)
  return PackedPoint2Point2Range(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Point2Point2Range}, d::PackedPoint2Point2Range)
  return Point2Point2Range(convert(SamplableBelief, d.Z))
end



"""
    $TYPEDEF

Range only measurement from Pose2 to Point2 variable.
"""
Base.@kwdef struct Pose2Point2Range{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T
  partial::Tuple{Int,Int} = (1,2)
end
Pose2Point2Range(Z::T) where {T <: IIF.SamplableBelief} = Pose2Point2Range(;Z)

getManifold(::Pose2Point2Range) = TranslationGroup(1)


function (cfo::CalcFactor{<:Pose2Point2Range})(rho, xi::Manifolds.ArrayPartition, lm)
  # Basically `EuclidDistance`
  return rho .- norm(lm .- xi.x[1])
end
# function (cfo::CalcFactor{<:Pose2Point2Range})(rho, xi::ProductRepr, lm)
#   # Basically `EuclidDistance`
#   return rho .- norm(lm .- xi.parts[1])
# end


Base.@kwdef struct PackedPose2Point2Range  <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{PackedPose2Point2Range}, d::Pose2Point2Range)
  return PackedPose2Point2Range(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Pose2Point2Range}, d::PackedPose2Point2Range)
  return Pose2Point2Range(convert(SamplableBelief, d.Z))
end
