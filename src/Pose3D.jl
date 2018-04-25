struct Pose3 <: IncrementalInference.InferenceVariable
  dims::Int
  Pose3() = new(6)
end
struct Point3 <: IncrementalInference.InferenceVariable
  dims::Int
  Point3() = new(3)
end
