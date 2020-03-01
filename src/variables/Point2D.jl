
"""
$(TYPEDEF)

XY Euclidean manifold variable node softtype.
"""
struct Point2 <: IncrementalInference.InferenceVariable
  dims::Int
  manifolds::Tuple{Symbol,Symbol}
  Point2() = new(2, (:Euclid, :Euclid))
end
