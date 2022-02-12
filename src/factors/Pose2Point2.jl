
export Pose2Point2, PackedPose2Point2

#-------------------------------------------------------------------------------
#

"""
    $TYPEDEF

Bearing and Range constraint from a Pose2 to Point2 variable.
"""
Base.@kwdef struct Pose2Point2{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
    Z::T = MvNormal(zeros(2),LinearAlgebra.diagm([0.01;0.01]))
    partial::Tuple{Int,Int} = (1,2)
end
# convenience and default constructor
Pose2Point2(x1::T=MvNormal(zeros(2),LinearAlgebra.diagm([0.01;0.01]))) where {T <: IIF.SamplableBelief} = Pose2Point2(Z=x1)

# TODO verify this is right for partial factor
getManifold(::InstanceType{Pose2Point2}) = getManifold(Point2)

# define the conditional probability constraint
function (cfo::CalcFactor{<:Pose2Point2})(meas,
                                          wXi,
                                          wLj )
  #
  #FIXME upgrade to manifolds 
  _wXi = getCoordinates(Pose2, wXi)
  wLj_pred = SE2(_wXi)*SE2([meas;0.0])
  return  wLj .- se2vee(wLj_pred)[1:2]
end


## Serialization support

Base.@kwdef struct PackedPose2Point2 <: AbstractPackedFactor
    Z::PackedSamplableBelief
end

function convert(::Type{PackedPose2Point2}, obj::Pose2Point2{T}) where {T <: IIF.SamplableBelief}
  return PackedPose2Point2(convert(PackedSamplableBelief, obj.Z))
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert(::Type{Pose2Point2}, packed::PackedPose2Point2)
  Pose2Point2(convert(SamplableBelief, packed.Z))
end
