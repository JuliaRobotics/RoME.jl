
"""
$(TYPEDEF)

XYZ Euclidean manifold variable node softtype.

Example
-------
```julia
p3 = Point3()
```
"""
struct Point3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol, Symbol, Symbol}
  Point3() = new(3, String[], (:Euclid,:Euclid,:Euclid))
end
