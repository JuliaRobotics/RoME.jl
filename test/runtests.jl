#addprocs(2)
#@show nprocs()

using RoME
using Test


@error("must restore testG2oParser.jl")
# "testG2oParser.jl";  ]

testfiles = [
"testTreeInitCommonMsg_IIF913.jl";
"testPose3Pose3NH.jl";
"testBeehive2D_CliqByCliq.jl";      # special case debugging
"testhigherdimroots.jl";
"testManifoldsPose2Equivalent.jl";
"testDidsonFunctions.jl";
"testPoint2Point2.jl";
"testPoint2Point2Init.jl";
"testBasicPose2Stationary.jl";
"TestPoseAndPoint2Constraints.jl";
"testDynPoint2D.jl";
"testBearingRange2D.jl";
"testDeltaOdo.jl";
"testFixedLagFG.jl";
"testMultimodalRangeBearing.jl";
"testDynPose2D.jl";
"threeDimLinearProductTest.jl";
"testPartialXYH.jl";
"testpartialpose3.jl";
"testpackingconverters.jl";
"TestDefaultFGInitialization.jl";
"testAccumulateFactors.jl";
"testDeadReckoningTether.jl"; 
"testPoint2Point2.jl";
"testParametric.jl";
"testInflation380.jl"]

## Tests not ready yet
# "HexagonalLightGraphs.jl"
# "testCameraFunctions.jl"
# "testmultiplefeatures.jl"
# "testPoint2Point2WorldBearing.jl";  # deprecate



for testf in testfiles
  println("[TEST] $testf")
  include(testf)
  println("[SUCCESS]")
  println()
  println()
  println()
end
