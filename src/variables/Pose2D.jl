

"""
$(TYPEDEF)

Pose2 is a SE(2) mechanization of two Euclidean translations and one Circular rotation, used for general 2D SLAM.
"""
struct Pose2 <: IncrementalInference.InferenceVariable
  dims::Int
  manifolds::Tuple{Symbol,Symbol,Symbol}
  Pose2() = new(3, (:Euclid, :Euclid, :Circular))
end
