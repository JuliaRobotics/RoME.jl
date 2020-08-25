"""
$(TYPEDEF)

XY Euclidean manifold variable node softtype.
"""
@defVariable Point2 2 (:Euclid, :Euclid)


"""
$(TYPEDEF)

XYZ Euclidean manifold variable node softtype.

Example
-------
```julia
p3 = Point3()
```
"""
@defVariable Point3 3 (:Euclid,:Euclid,:Euclid)


"""
$(TYPEDEF)

Pose2 is a SE(2) mechanization of two Euclidean translations and one Circular rotation, used for general 2D SLAM.
"""
@defVariable Pose2 3 (:Euclid,:Euclid,:Circular)

"""
$(TYPEDEF)

Pose3 is currently a Euler angle mechanization of three Euclidean translations and three Circular rotation.

Future:
------
- Work in progress on AMP3D for proper non-Euler angle on-manifold operations.
- TODO the AMP upgrade is aimed at resolving 3D to Quat/SE3/SP3 -- current Euler angles will be replaced
"""
@defVariable Pose3 6 (:Euclid,:Euclid,:Euclid,:Circular,:Circular,:Circular)

"""
$(TYPEDEF)

Dynamic point in 2D space with velocity components: `x, y, dx/dt, dy/dt`

"""
@defVariable DynPoint2 4 (:Euclid,:Euclid,:Euclid,:Euclid)

"""
$(TYPEDEF)

Dynamic pose variable with velocity components: `x, y, theta, dx/dt, dy/dt`
"""
@defVariable DynPose2 5 (:Euclid,:Euclid,:Circular,:Euclid,:Euclid)