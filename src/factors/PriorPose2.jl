
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
Base.@kwdef struct PriorPose2{T <: SamplableBelief} <: IIF.AbstractPrior
  Z::T = MvNormal(zeros(3), diagm([1;1;0.1]))
end

DFG.getManifold(::InstanceType{PriorPose2}) = getManifold(Pose2) # SpecialEuclidean(2)

@inline function _vee(::SpecialEuclidean{2}, X::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}) where T<:Real
  return SVector{3,T}(X.x[1][1],X.x[1][2],X.x[2][2])
end

@inline function _compose(::SpecialEuclidean{2}, p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, q::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}) where T<:Real
  return ArrayPartition(p.x[1] + p.x[2]*q.x[1], p.x[2]*q.x[2])
end

function (cf::CalcFactor{<:PriorPose2})(_m::AbstractArray{MT}, _p::AbstractArray{PT})  where {MT<:Real,PT<:Real}
  T = promote_type(MT, PT)
  m = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _m)
  p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
  return cf(m,p)
end

# TODO the log here looks wrong (for gradients), consider:
# X = log(p⁻¹ ∘ m) 
# X = log(M, ϵ, Manifolds.compose(M, inv(M, p), m))
function (cf::CalcFactor{<:PriorPose2})(
            m::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, 
            p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}) where T<:Real

  M = getManifold(Pose2)
  Xc = _vee(M, log(M, p, m))
  return Xc
end

#TODO Serialization of reference point p 
## Serialization support

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPriorPose2  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  return PackedPriorPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  return PriorPose2(convert(SamplableBelief, d.Z))
end




## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  compareDensity(a.Z, b.Z)
end
