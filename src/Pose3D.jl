"""
$(TYPEDEF)
"""
struct Point3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol, Symbol, Symbol}
  Point3() = new(3, String[], (:Euclid,:Euclid,:Euclid))
end

"""
$(TYPEDEF)
"""
struct Pose3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  manifolds::Tuple{Symbol,Symbol,Symbol,Symbol,Symbol,Symbol}
  # TODO the AMP upgrade is aimed at resolving 3D to Quat/SE3/SP3 -- current Euler angles will be replaced
  Pose3() = new(6, String[], (:Euclid,:Euclid,:Euclid,:Euclid,:Euclid,:Euclid) )
end
