
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

DFG.getManifold(::InstanceType{Pose2Pose2}) = LeftInvariantMetricSE(2)

Pose2Pose2(::UniformScaling) = Pose2Pose2()

function (cf::CalcFactor{<:Pose2Pose2})(X, p, q)
    # X ∈ TₚM, X̂ ∈ TₚM, p,q ∈ M
    M = getManifold(Pose2Pose2)
    X̂ = log(M, p, q)
    return vee(M, p, X - X̂) # TODO check sign
end

# An alternative error function that parallel transports the error vector back to p
# function (cf::CalcFactor{<:RobustPose2Pose2})(X, p, q)
#     G = getManifold(RobustPose2Pose2)
#     q̂ = exp(G, p, X)
#     E_q̂ = log(G, q̂, q)
#     # Parallel transport the error vector back to p.
#     E_p = parallel_transport_to(G, q̂, E_q̂, p)    
#     return vee(LieAlgebra(G), E_p) 
# end


# NOTE, serialization support -- will be reduced to macro in future
# ------------------------------------

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPose2Pose2 <: AbstractPackedObservation
    Z::PackedBelief
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
    return Pose2Pose2(convert(SamplableBelief, d.Z))
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
    return PackedPose2Pose2(convert(PackedBelief, d.Z))
end

# FIXME, rather have separate compareDensity functions
function compare(a::Pose2Pose2, b::Pose2Pose2; tol::Float64 = 1e-10)
    return compareDensity(a.Z, b.Z)
end

#
