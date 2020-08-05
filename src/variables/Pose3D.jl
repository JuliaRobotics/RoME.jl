
"""
$(TYPEDEF)

Pose3 is currently a Euler angle mechanization of three Euclidean translations and three Circular rotation.

Future:
------
- Work in progress on AMP3D for proper non-Euler angle on-manifold operations.
"""
struct Pose3 <: IncrementalInference.InferenceVariable end
getDimension(::Pose3) = 6
# TODO the AMP upgrade is aimed at resolving 3D to Quat/SE3/SP3 -- current Euler angles will be replaced
getManifolds(::Pose3) = (:Euclid,:Euclid,:Euclid,:Circular,:Circular,:Circular)
