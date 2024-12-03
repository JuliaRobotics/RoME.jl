
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
Pose2Point2(Z::SamplableBelief) = Pose2Point2(;Z)

# TODO verify this is right for partial factor
getManifold(::InstanceType{Pose2Point2}) = getManifold(Point2)

# define the conditional probability constraint
function (cfo::CalcFactor{<:Pose2Point2})(p_Xpq,
                                          w_T_p,
                                          w_Tl_q )
  #
  M = SpecialEuclidean(2; vectors=HybridTangentRepresentation())

  p_T_qhat = ArrayPartition(SA[p_Xpq[1];p_Xpq[2]], SMatrix{2,2}([1 0; 0 1.]))
  _w_T_p = ArrayPartition(SA[w_T_p.x[1]...], SMatrix{2,2}(w_T_p.x[2]))
  w_H_qhat = affine_matrix(M, _w_T_p) * affine_matrix(M, p_T_qhat)
  # Issue, this produces ComposedFunction which errors on not having field .x[1]
  # w_T_qhat = compose(M, _w_T_p, p_T_qhat)
  return w_Tl_q .- w_H_qhat[1:2,end]
  
  # upgraded to manifolds 
  # w_Cp = getCoordinates(Pose2, w_T_p)
  # w_Tl_q_pred = SE2(w_Cp)*SE2([p_Xpq;0.0])
  # return  w_Tl_q .- se2vee(w_Tl_q_pred)[1:2]
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
