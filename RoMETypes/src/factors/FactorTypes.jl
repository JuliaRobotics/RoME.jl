# ========================================================================================
# Pose
# ========================================================================================

"""
$(TYPEDEF)

Direction observation information of a `Point3` variable.
"""
DFG.@defObservationType PriorPoint3 PriorObservation LieGroups.TranslationGroup(3)
DFG.@defObservationType Point3Point3 RelativeObservation LieGroups.TranslationGroup(3)


DFG.@defObservationType Point2Point2 RelativeObservation LieGroups.TranslationGroup(2)
DFG.@defObservationType PriorPoint2 PriorObservation LieGroups.TranslationGroup(2)

# ========================================================================================
# Pose
# ========================================================================================

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
DFG.@defObservationType Pose2Pose2 RelativeObservation SOnxRn_MetricManifold(2)

"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
DFG.@defObservationType PriorPose2 PriorObservation TranslationGroup(2) × SpecialOrthogonalGroup(2)

DFG.@defObservationType Pose3Pose3 RelativeObservation SOnxRn_MetricManifold(3)
DFG.@defObservationType PriorPose3 PriorObservation TranslationGroup(3) × SpecialOrthogonalGroup(3)
