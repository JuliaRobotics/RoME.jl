
"""
$(TYPEDEF)

Rigid transform between two Pose2's, assuming (x,y,theta).

Calcuated as:
```math
\\begin{aligned}
\\hat{X}_p=\\log_pq\\\\
X^i = \\mathrm{vee}(X_m-\\hat{X})
\\end{aligned}
```
with:
``\\mathcal M= \\mathrm{SE}(2)`` Special Euclidean group\\
``p`` and ``q`` ``\\in \\mathcal M`` the two Pose2 points\\
the measurement vector ``X_m \\in T_p \\mathcal M``\\
and the error vector ``\\hat{X} \\in T_p \\mathcal M``\\
``X^i`` coordinates of ``X``

DevNotes
- Maybe with Manifolds.jl, `{T <: IIF.SamplableBelief, S, R, P}`

Related

[`Pose3Pose3`](@ref), [`Point2Point2`](@ref), [`MutablePose2Pose2Gaussian`](@ref), [`DynPose2`](@ref), [`IMUDeltaFactor`](@ref)
"""
Base.@kwdef struct Pose2Pose2{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
    Z::T = MvNormal(Diagonal([1.0; 1.0; 1.0]))
end

DFG.getManifold(::InstanceType{Pose2Pose2}) = SOnxRn_MetricManifold(2)

Pose2Pose2(::UniformScaling) = Pose2Pose2()

function (cf::CalcFactor{<:Pose2Pose2})(X, p, q)
    # X ∈ TₚM, X̂ ∈ TₚM, p,q ∈ M
    M = getManifold(Pose2Pose2)
    X̂ = log(M, p, q)
    return vee(M, p, X - X̂) # TODO check sign
end

# NOTE, serialization support -- will be reduced to macro in future
# ------------------------------------

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPose2Pose2 <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
    return Pose2Pose2(convert(SamplableBelief, d.Z))
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
    return PackedPose2Pose2(convert(PackedSamplableBelief, d.Z))
end

# FIXME, rather have separate compareDensity functions
function compare(a::Pose2Pose2, b::Pose2Pose2; tol::Float64 = 1e-10)
    return compareDensity(a.Z, b.Z)
end

#
