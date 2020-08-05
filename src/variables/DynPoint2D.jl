
"""
$(TYPEDEF)

Dynamic point in 2D space with velocity components: `x, y, dx/dt, dy/dt`

"""
mutable struct DynPoint2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  manifolds::Tuple{Symbol, Symbol, Symbol, Symbol}
  DynPoint2(;ut::Int64=-9999999999) = new(ut, 4,(:Euclid,:Euclid,:Euclid,:Euclid,))
end
getDimension(::DynPoint2) = 4
getManifolds(::DynPoint2) = (:Euclid,:Euclid,:Euclid,:Euclid)