
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
  manifolds::Tuple{Symbol, Symbol, Symbol}
  Point3() = new(3, (:Euclid,:Euclid,:Euclid))
end
