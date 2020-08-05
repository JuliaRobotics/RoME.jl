
"""
$(TYPEDEF)

XY Euclidean manifold variable node softtype.
"""
struct Point2 <: IncrementalInference.InferenceVariable end
getDimension(::Point2) = 2
getManifolds(::Point2) = (:Euclid, :Euclid)

