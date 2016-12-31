addprocs(3)
@show nprocs()

using RoME
using Base.Test

using JLD, HDF5



println("[TEST] Camera function evaluations...")
include("testCameraFunctions.jl")
println("[SUCCESS]")

println("[TEST] Linear array functions evaluations...")
include("testDidsonFunctions.jl")
println("[SUCCESS]")

println("[TEST] Pose2 evaluations...")
include("TestPoseAndPoint2Constraints.jl")
println("[SUCCESS]")

println("[TEST] Pose3 evaluations...")
include("threeDimLinearProductTest.jl")
println("[SUCCESS]")



println("[TEST] ensure Pose3Pose3NH evaluations...")
include("testPose3Pose3NH.jl")
println("[SUCCESS]")
