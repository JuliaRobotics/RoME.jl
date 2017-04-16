addprocs(2)
@show nprocs()

using RoME
using Base.Test

using JLD, HDF5

println("[TEST] numeric root...")
include("testhigherdimroots.jl")
println("[SUCCESS]")

println("[TEST] Camera function evaluations...")
include("testCameraFunctions.jl")
println("[SUCCESS]")

println("[TEST] MultipleFeatures constraints")
include("testmultiplefeatures.jl")
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

println("[TEST] saving to and loading from .jld file")
savejld(fg, file="tempfg.jld" )
fgu = loadjld( file="tempfg.jld" )
Base.rm("tempfg.jld")
println("Success")

println("[TEST] partial pose3 evaluations...")
include("testpartialpose3.jl")
println("[SUCCESS]")

println("[TEST] PartialPose3XYYaw evaluations...")
include("testPartialXYH.jl")
println("[SUCCESS]")

println("[TEST] packing converters...")
include("testpackingconverters.jl")
println("[SUCCESS]")
