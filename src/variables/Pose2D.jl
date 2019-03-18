

"""
$(TYPEDEF)

Pose2 is a SE(2) mechanization of two Euclidean translations and one Circular rotation, used for general 2D SLAM.
"""
struct Pose2 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol,Symbol,Symbol}
  Pose2(;labels::Vector{<:AbstractString}=String[]) = new(3, labels, (:Euclid, :Euclid, :Circular))
end
