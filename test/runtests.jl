#addprocs(2)
#@show nprocs()

using RoME
using Test

# might be unnecessary
# using JLD2  #, HDF5

@error "not testing "*"testPoint2Point2WorldBearing.jl";

testfiles = [ "testhigherdimroots.jl";
"testManifoldsPose2Equivalent.jl";
"testDidsonFunctions.jl";
"testPoint2Point2.jl";
"testPoint2Point2Init.jl";
"testBasicPose2Stationary.jl";
"TestPoseAndPoint2Constraints.jl";
"testBeehive2D_CliqByCliq.jl";
"testDynPoint2D.jl";
"testBearingRange2D.jl";
"testDeltaOdo.jl";
"testFixedLagFG.jl";
"testMultimodalRangeBearing.jl";
"testDynPose2D.jl";
"threeDimLinearProductTest.jl";
"testPose3Pose3NH.jl";
"testPartialXYH.jl";
"testpartialpose3.jl";
"testpackingconverters.jl";
"TestDefaultFGInitialization.jl";
"testAccumulateFactors.jl";
"testDeadReckoningTether.jl";
"testG2oParser.jl";  ]
# "HexagonalLightGraphs.jl"
# "testCameraFunctions.jl"
# "testmultiplefeatures.jl"


for testf in testfiles
  println("[TEST] $testf")
  include(testf)
  println("[SUCCESS]")
  println()
  println()
  println()
end

    #
    # println("[TEST] numeric root...")
    # include("testhigherdimroots.jl")
    # println("[SUCCESS]")
    #
    # include("testManifoldsPose2Equivalent.jl")
    #
    # include("testDidsonFunctions.jl")
    #
    # include("testPoint2Point2.jl")
    #
    # include("testPoint2Point2Init.jl")
    #
    # include("testBasicPose2Stationary.jl")
    #
    # include("TestPoseAndPoint2Constraints.jl")
    #
    # include("testBeehive2D_CliqByCliq.jl")
    #
    # include("testDynPoint2D.jl")
    #
    # include("testPoint2Point2WorldBearing.jl")
    #
    # include("testBearingRange2D.jl")
    #
    # include("testDeltaOdo.jl")
    #
    # include("testFixedLagFG.jl")
    #
    # include("testMultimodalRangeBearing.jl")
    #
    # include("testDynPose2D.jl")
    #
    # println("[TEST] Pose3 evaluations...")
    # include("threeDimLinearProductTest.jl")
    # println("[SUCCESS]")
    #
    # println("[TEST] ensure Pose3Pose3NH evaluations...")
    # include("testPose3Pose3NH.jl")
    # println("[SUCCESS]")
    #
    # println("[TEST] PartialPose3XYYaw evaluations...")
    # include("testPartialXYH.jl")
    # println("[SUCCESS]")
    #
    # println("[TEST] partial pose3 evaluations...")
    # include("testpartialpose3.jl")
    # println("[SUCCESS]")
    #
    # println("[TEST] packing converters...")
    # include("testpackingconverters.jl")
    # println("[SUCCESS]")
    #
    # include("TestDefaultFGInitialization.jl")
    #
    # include("testAccumulateFactors.jl")
    #
    # include("testDeadReckoningTether.jl")
    #
    # include("testG2oParser.jl")

## not ready for use
## include("HexagonalLightGraphs.jl")

# Requires standardized testing strategy
# println("[TEST] Camera function evaluations...")
# include("testCameraFunctions.jl")
# println("[SUCCESS]")

@warn "Skipping multiple feature constraint test for the time being"
# println("[TEST] MultipleFeatures constraints")
# include("testmultiplefeatures.jl")
# println("[SUCCESS]")
