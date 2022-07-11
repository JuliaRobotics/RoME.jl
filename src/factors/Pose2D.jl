

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

[`Pose3Pose3`](@ref), [`Point2Point2`](@ref), [`MutablePose2Pose2Gaussian`](@ref), [`DynPose2`](@ref), [`InertialPose3`](ref)
"""
Base.@kwdef struct Pose2Pose2{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T = MvNormal(Diagonal([1.0; 1.0; 1.0]))
end

DFG.getManifold(::InstanceType{Pose2Pose2}) = getManifold(Pose2) # Manifolds.SpecialEuclidean(2)

Pose2Pose2(::UniformScaling) = Pose2Pose2()

function preambleCache(dfg::AbstractDFG, vars::AbstractVector{<:DFGVariable}, pp::Pose2Pose2)
  M = getManifold(pp)
  (;manifold=M, ϵ0=getPointIdentity(M), Xc=zeros(3), q̂=getPointIdentity(M))
end

# Assumes X is a tangent vector
function (cf::CalcFactor{<:Pose2Pose2})(X, p, q)

    q̂ = allocate(q) 
    M = getManifold(Pose2)
    ϵ0 = getPointIdentity(M)
    exp!(M, q̂, ϵ0, X)
    Manifolds.compose!(M, q̂, p, q̂)   
    Xc = vee(M, q, log!(M, q̂, q, q̂))
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
