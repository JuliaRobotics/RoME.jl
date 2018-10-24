#addprocs(2)
#@show nprocs()

using RoME
using Test

# might be unnecessary
using JLD2  #, HDF5

println("[TEST] numeric root...")
include("testhigherdimroots.jl")
println("[SUCCESS]")

# Requires standardized testing strategy
# println("[TEST] Camera function evaluations...")
# include("testCameraFunctions.jl")
# println("[SUCCESS]")

println("[TEST] MultipleFeatures constraints")
include("testmultiplefeatures.jl")
println("[SUCCESS]")

println("[TEST] Linear array function evaluations...")
include("testDidsonFunctions.jl")
println("[SUCCESS]")

include("testPoint2Point2.jl")

println("[TEST] Pose2 evaluations...")
include("TestPoseAndPoint2Constraints.jl")
println("[SUCCESS]")


@testset "[TEST] Pose2 evaluations..." begin
  include("testDynPoint2D.jl")
end

include("testPoint2Point2WorldBearing.jl")

include("testBearingRange2D.jl")

include("testFixedLagFG.jl")

include("testMultimodalRangeBearing.jl")

include("testDynPose2D.jl")

println("[TEST] Pose3 evaluations...")
include("threeDimLinearProductTest.jl")
println("[SUCCESS]")

println("[TEST] ensure Pose3Pose3NH evaluations...")
include("testPose3Pose3NH.jl")
println("[SUCCESS]")

println("[TEST] saving to and loading from .jld2 file")
savejld(fg, file="tempfg.jld2" )
fgu = loadjld( file="tempfg.jld2" )
Base.rm("tempfg.jld2")
println("Success")

println("[TEST] PartialPose3XYYaw evaluations...")
include("testPartialXYH.jl")
println("[SUCCESS]")

println("[TEST] partial pose3 evaluations...")
include("testpartialpose3.jl")
println("[SUCCESS]")

println("[TEST] packing converters...")
include("testpackingconverters.jl")
println("[SUCCESS]")

include("TestDefaultFGInitialization.jl")
