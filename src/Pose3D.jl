struct Pose3 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  Pose3() = new(6, String["POSE";])
end
struct Point3 <: IncrementalInference.InferenceVariable
  dims::Int
  Point3() = new(3)
end
