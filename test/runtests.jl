#addprocs(2)
#@show nprocs()

using RoME
using Test

@testset "Numeric Root" begin
include("testhigherdimroots.jl")
end

@testset "Manifolds" begin
include("testManifoldsPose2Equivalent.jl")
end

# Requires standardized testing strategy
# println("[TEST] Camera function evaluations...")
# include("testCameraFunctions.jl")
# end

@testset "MultipleFeatures constraints" begin
include("testmultiplefeatures.jl")
end

@testset "Linear array function evaluations" begin
include("testDidsonFunctions.jl")
end

@testset "Point2Point2 Tests" begin
include("testPoint2Point2.jl")
end

@testset "Point2Point2 Init Tests" begin
include("testPoint2Point2Init.jl")
end

@testset "BasicPose2Stationary Tests" begin
include("testBasicPose2Stationary.jl")
end

@testset "Pose2 Evaluations" begin
include("TestPoseAndPoint2Constraints.jl")
end

@testset "DFG End-to-End Tests" begin
include("dfg/HexagonalGraphs.jl")
include("dfg/FileDFG.jl")
include("dfg/HexagonalLightGraphs.jl")
# Don't run the cloud tests yet
end

@testset "Beehive Tests" begin
include("testBeehive2D_CliqByCliq.jl")
end

@testset "Pose2 evaluations..." begin
include("testDynPoint2D.jl")
end

@testset "Point2Point2WorldBearing Tests" begin
include("testPoint2Point2WorldBearing.jl")
end

@testset "BearingRange2D Tests" begin
include("testBearingRange2D.jl")
end

@testset "FixedLag Tests" begin
include("testFixedLagFG.jl")
end

@testset "MultimodalRangeBearing Tests" begin
include("testMultimodalRangeBearing.jl")
end
@testset "DynPose2D Tests" begin
include("testDynPose2D.jl")
end

@testset "Pose3 evaluations" begin
include("threeDimLinearProductTest.jl")
end

@testset "Ensure Pose3Pose3NH evaluations" begin
include("testPose3Pose3NH.jl")
end

@testset "PartialPose3XYYaw evaluations" begin
include("testPartialXYH.jl")
end

@testset "Partial pose3 evaluations" begin
include("testpartialpose3.jl")
end

@testset "Packing converters" begin
include("testpackingconverters.jl")
end

@testset "Default FG Initialization" begin
include("TestDefaultFGInitialization.jl")
end
