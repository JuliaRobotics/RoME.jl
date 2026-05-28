"""
$(TYPEDEF)

XY Euclidean manifold variable node softtype.
"""
@defStateType Point2 TranslationGroup(2) SA[0.0; 0.0]

"""
$(TYPEDEF)

XYZ Euclidean manifold variable node softtype.

Example
-------
```julia
p3 = Point3()
```
"""
@defStateType Point3 TranslationGroup(3) SA[0; 0; 0.0]

"""
$(TYPEDEF)

Pose2 is a SE(2) mechanization of two Euclidean translations and one Circular rotation, used for general 2D SLAM.
"""
@defStateType(
    Pose2,
    TranslationGroup(2) × SpecialOrthogonalGroup(2), #TODO look at using LeftInvariantMetricSE(2) 
    ArrayPartition(SA[0; 0.0], SA[1 0; 0 1.0])
)

"""
$(TYPEDEF)

Pose3 is currently a Euler angle mechanization of three Euclidean translations and three Circular rotation.

Future:
------
- Work in progress on AMP3D for proper non-Euler angle on-manifold operations.
- TODO the AMP upgrade is aimed at resolving 3D to Quat/SE3/SP3 -- current Euler angles will be replaced
"""
@defStateType(
    Pose3,
    TranslationGroup(3) × SpecialOrthogonalGroup(3), #TODO look at using LeftInvariantMetricSE(3) 
    ArrayPartition(SA[0; 0; 0.0], SA[1 0 0; 0 1 0; 0 0 1.0])
)

@defStateType Rotation3 SpecialOrthogonalGroup(3) SA[1 0 0; 0 1 0; 0 0 1.0]

@defStateType(
    RotVelPos,
    SpecialOrthogonalGroup(3) × TranslationGroup(3) × TranslationGroup(3),
    ArrayPartition(SA[1 0 0; 0 1 0; 0 0 1.0], SA[0; 0; 0.0], SA[0; 0; 0.0])
)

# 3 translations and 3 velocity in graph-base-frame
@defStateType(
    VelPos3,
    TranslationGroup(3) × TranslationGroup(3),
    ArrayPartition(SA[0; 0; 0.0], SA[0; 0; 0.0])
)

# @defStateType VelPose3 Manifolds.ProductGroup(ProductManifold(TranslationGroup(3), TranslationGroup(3), SpecialOrthogonal(3))) ArrayPartition(SA[0; 0; 0.0], SA[0;0;0.0], SA[1 0 0; 0 1 0; 0 0 1.0])
# Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{typeof(getManifold(VelPose3))}) =
#     (:Euclid,:Euclid,:Euclid,:Euclid,:Euclid,:Euclid, :Circular,:Circular,:Circular,)

"""
$(TYPEDEF)

Dynamic point in 2D space with velocity components: `x, y, dx/dt, dy/dt`

"""
@defStateType DynPoint2 TranslationGroup(4) zero(SVector{4, Float64})

"""
$(TYPEDEF)

Dynamic pose variable with velocity components: `x, y, theta, dx/dt, dy/dt`

Note
- The `SE2E2_Manifold` definition used currently is a hack to simplify the transition to Manifolds.jl, see #244 
- Replaced `SE2E2_Manifold` hack with `ProductManifold(SpecialEuclidean(2), TranslationGroup(2))`, confirm if it is correct.
"""
@defStateType(
    DynPose2,
    # LeftInvariantMetricSE(2) × TranslationGroup(2), #FIXME SOnxRn(2) or SE(2)
    # ArrayPartition(ArrayPartition(SA[0;0.0],SA[1 0; 0 1.0]),SA[0;0.0])
    TranslationGroup(2) × SpecialOrthogonalGroup(2) × TranslationGroup(2),
    ArrayPartition(SA[0; 0.0], SA[1 0; 0 1.0], SA[0; 0.0])
)

"""
$SIGNATURES

Function to project only XY data onto Cartesian plane for 2D plotting.
"""
projectCartesian(
    pose::Union{<:Point2, <:Point3, <:Pose2, <:Pose3, <:DynPoint2, <:DynPose2},
    x::Vector{Float64},
) = [x[1]; x[2]; 0]

#
