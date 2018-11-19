"""
$(TYPEDEF)
"""
struct Point3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  Point3() = new(3, String[])
end

"""
$(TYPEDEF)
"""
struct Pose3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  Pose3() = new(6, String["POSE";])
end
