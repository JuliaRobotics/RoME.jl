# TODO integration underway with Manifolds.jl, see RoME #244, also see IIF #467 regarding consolidation effort.

export projectCartesian


"""
$(TYPEDEF)

XY Euclidean manifold variable node softtype.
"""
@defVariable Point2 Euclidean(2)


"""
$(TYPEDEF)

XYZ Euclidean manifold variable node softtype.

Example
-------
```julia
p3 = Point3()
```
"""
@defVariable Point3 Euclidean(3)


"""
$(TYPEDEF)

Pose2 is a SE(2) mechanization of two Euclidean translations and one Circular rotation, used for general 2D SLAM.
"""
@defVariable Pose2 SpecialEuclidean(2)

"""
$(TYPEDEF)

Pose3 is currently a Euler angle mechanization of three Euclidean translations and three Circular rotation.

Future:
------
- Work in progress on AMP3D for proper non-Euler angle on-manifold operations.
- TODO the AMP upgrade is aimed at resolving 3D to Quat/SE3/SP3 -- current Euler angles will be replaced
"""
@defVariable Pose3 SpecialEuclidean(3)

"""
$(TYPEDEF)

Dynamic point in 2D space with velocity components: `x, y, dx/dt, dy/dt`

"""
@defVariable DynPoint2 Euclidean(4)

"""
$(TYPEDEF)

Dynamic pose variable with velocity components: `x, y, theta, dx/dt, dy/dt`

Note
- The `SE2E2_Manifold` definition used currently is a hack to simplify the transition to Manifolds.jl, see #244 
"""
@defVariable DynPose2 SE2E2_Manifold



"""
$SIGNATURES

Function to project only XY data onto Cartesian plane for 2D plotting.
"""
projectCartesian(pose::Union{<:Point2,<:Point3,<:Pose2,<:Pose3,<:DynPoint2,<:DynPose2}, 
                 x::Vector{Float64}) = [x[1]; x[2]; 0]

#



# Still experimental
# export BearingRange2
@defVariable BearingRange2 BearingRange_Manifold




#
