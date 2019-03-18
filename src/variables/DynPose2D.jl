"""
$(TYPEDEF)

Dynamic pose variable with velocity components: `x, y, theta, dx/dt, dy/dt`
"""
mutable struct DynPose2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol,Symbol,Symbol,Symbol,Symbol}
  DynPose2(;ut::Int64=-9999999999, labels::Vector{<:AbstractString}=String[]) = new(ut, 5, labels, (:Euclid,:Euclid,:Circular,:Euclid,:Euclid))
end
