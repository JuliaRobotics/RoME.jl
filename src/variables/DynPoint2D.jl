
"""
$(TYPEDEF)

Dynamic point in 2D space with velocity components: `x, y, dx/dt, dy/dt`

"""
mutable struct DynPoint2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol, Symbol, Symbol, Symbol}
  DynPoint2(;ut::Int64=-9999999999, labels::Vector{<:AbstractString}=String[]) = new(ut, 4, labels,(:Euclid,:Euclid,:Euclid,:Euclid,))
end
