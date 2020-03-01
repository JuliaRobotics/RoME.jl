"""
$(TYPEDEF)

Dynamic pose variable with velocity components: `x, y, theta, dx/dt, dy/dt`
"""
mutable struct DynPose2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  manifolds::Tuple{Symbol,Symbol,Symbol,Symbol,Symbol}
  DynPose2(;ut::Int64=-9999999999) = new(ut, 5, (:Euclid,:Euclid,:Circular,:Euclid,:Euclid))
end
