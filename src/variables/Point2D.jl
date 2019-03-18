
"""
$(TYPEDEF)

XY Euclidean manifold variable node softtype.
"""
struct Point2 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol,Symbol}
  Point2(;labels::Vector{<:AbstractString}=String[]) = new(2, labels, (:Euclid, :Euclid))
end
