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

Rigid transform between two Pose2's.

Calculated as:
```math
\\begin{aligned}
\\hat{X} = \\log(\\mathcal{M}, p, q)\\\\
X_e = X_m - \\hat{X}\\\\
X^i = \\mathrm{vee}(\\mathfrak{g}, X_e)
\\end{aligned}
```
with:
`\\mathcal M= LeftInvariantMetricSE(2)` Special Euclidean group with a left-invariant metric\
`\\mathfrak{g} = \\mathfrak{se}(2)` the Lie algebra at the identity element\
`p` and `q` `\\in \\mathcal M` the two Pose2 points\
the measurement vector `X_m \\in \\mathfrak{g}`\
the predicted relative tangent vector `\\hat{X} \\in \\mathfrak{g}`\
the error vector `X_e \\in \\mathfrak{g}`\
`X^i` coordinate vector of the error

Related

[`Pose3Pose3`](@ref), [`Point2Point2`](@ref), [`MutablePose2Pose2Gaussian`](@ref), [`DynPose2`](@ref), [`IMUDeltaFactor`](@ref)
"""
DFG.@defObservationType Pose2Pose2 RelativeObservation LeftInvariantMetricSE(2)

"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
DFG.@defObservationType(
    PriorPose2,
    PriorObservation,
    TranslationGroup(2) × SpecialOrthogonalGroup(2)
)

DFG.@defObservationType Pose3Pose3 RelativeObservation LeftInvariantMetricSE(3)

DFG.@defObservationType(
    PriorPose3,
    PriorObservation,
    TranslationGroup(3) × SpecialOrthogonalGroup(3)
)
