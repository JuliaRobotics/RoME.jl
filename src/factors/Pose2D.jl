

"""
$(TYPEDEF)

Rigid transform between two Pose2's, assuming (x,y,theta).

Calcuated as:
```math
\\begin{aligned}
\\hat{q}=\\exp_pX_m\\\\
X = \\log_q \\hat{q}\\\\
X^i = \\mathrm{vee}(q, X)
\\end{aligned}
```
with:
``\\mathcal M= \\mathrm{SE}(2)`` Special Euclidean group\\
``p`` and ``q`` ``\\in \\mathcal M`` the two Pose2 points\\
the measurement vector ``X_m \\in T_p \\mathcal M``\\
and the error vector ``X \\in T_q \\mathcal M``\\
``X^i`` coordinates of ``X``

DevNotes
- Maybe with Manifolds.jl, `{T <: IIF.SamplableBelief, S, R, P}`

Related

[`Pose3Pose3`](@ref), [`Point2Point2`](@ref), [`MutablePose2Pose2Gaussian`](@ref), [`DynPose2`](@ref), [`IMUDeltaFactor`](@ref)
"""
Base.@kwdef struct Pose2Pose2{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T = MvNormal(Diagonal([1.0; 1.0; 1.0]))
end

DFG.getManifold(::InstanceType{Pose2Pose2}) = Manifolds.SpecialEuclidean(2)

Pose2Pose2(::UniformScaling) = Pose2Pose2()


# Assumes X is a tangent vector
function (cf::CalcFactor{<:Pose2Pose2})(_X::AbstractArray{MT}, _p::AbstractArray{PT}, _q::AbstractArray{LT})  where {MT,PT,LT}
  #TODO remove this convertions
  # @warn "This warning should not be triggered after StaticArrays upgrade" maxlog=10
  T = promote_type(MT, PT, LT)
  X = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _X)
  p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
  q = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _q)
  return cf(X,p,q)
end

# function calcPose2Pose2(
function (cf::CalcFactor{<:Pose2Pose2})(
              X::ArrayPartition{XT, Tuple{SVector{2, XT}, SMatrix{2, 2, XT, 4}}},
              p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, 
              q::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}) where {XT<:Real,T<:Real}

    M = getManifold(Pose2)
    # ϵ0 = ArrayPartition(zeros(SVector{2,T}), SMatrix{2, 2, T}(I))
    ϵ0 = getPointIdentity(M)

    ϵX = exp(M, ϵ0, X)
    # q̂ = Manifolds.compose(M, p, ϵX)    
    q̂ = _compose(M, p, ϵX)
    X_hat = log(M, q, q̂)#::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}
    # Xc = vee(M, q, X_hat)
    Xc = _vee(M, X_hat)#::SVector{3,T}
    return Xc
end


# NOTE, serialization support -- will be reduced to macro in future
# ------------------------------------

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPose2Pose2  <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
  return Pose2Pose2(convert(SamplableBelief, d.Z))
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
  return PackedPose2Pose2(convert(PackedSamplableBelief, d.Z))
end


# FIXME, rather have separate compareDensity functions
function compare(a::Pose2Pose2,b::Pose2Pose2; tol::Float64=1e-10)
  return compareDensity(a.Z, b.Z)
end

#
