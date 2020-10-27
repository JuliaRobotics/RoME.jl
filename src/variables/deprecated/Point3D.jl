
"""
$(TYPEDEF)

XYZ Euclidean manifold variable node softtype.

Example
-------
```julia
p3 = Point3()
```
"""
struct Point3 <: IncrementalInference.InferenceVariable end
getDimension(::Point3) = 3
getManifolds(::Point3) = (:Euclid,:Euclid,:Euclid)

