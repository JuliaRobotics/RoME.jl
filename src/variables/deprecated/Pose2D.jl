

"""
$(TYPEDEF)

Pose2 is a SE(2) mechanization of two Euclidean translations and one Circular rotation, used for general 2D SLAM.
"""
struct Pose2 <: IncrementalInference.InferenceVariable end
getDimension(::Pose2) = 3
getManifolds(::Pose2) = (:Euclid,:Euclid,:Circular)
