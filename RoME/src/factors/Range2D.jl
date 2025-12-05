
"""
$(TYPEDEF)
"""
DFG.@defObservationType Point2Point2Range RelativeObservation TranslationGroup(1)

function (cfo::CalcFactor{<:Point2Point2Range})(rho, xi, lm)
    # Basically `EuclidDistance`
    # must return all dimensions
    return rho .- norm(lm[1:2] .- xi[1:2])
end

"""
    $TYPEDEF

Range only measurement from Pose2 to Point2 variable.
"""
DFG.@kwarg struct Pose2Point2Range{T} <: IIF.AbstractManifoldMinimize
    Z::T & DFG.@packed
    partial::Tuple{Int, Int} = (1, 2)
end
Pose2Point2Range(Z::T) where {T <: IIF.SamplableBelief} = Pose2Point2Range(; Z)

DFG.getManifold(::Type{<:Pose2Point2Range}) = TranslationGroup(1)

function (cfo::CalcFactor{<:Pose2Point2Range})(rho, xi::ArrayPartition, lm)
    # Basically `EuclidDistance`
    return rho .- norm(lm .- xi.x[1])
end
