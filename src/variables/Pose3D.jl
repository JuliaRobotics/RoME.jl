
"""
$(TYPEDEF)

Pose3 is currently a Euler angle mechanization of three Euclidean translations and three Circular rotation.

Future:
------
- Work in progress on AMP3D for proper non-Euler angle on-manifold operations.
"""
struct Pose3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol,Symbol,Symbol,Symbol,Symbol,Symbol}
  # TODO the AMP upgrade is aimed at resolving 3D to Quat/SE3/SP3 -- current Euler angles will be replaced
  Pose3() = new(6, String[], (:Euclid,:Euclid,:Euclid,:Circular,:Circular,:Circular) )
end
