#addprocs(2)
#@show nprocs()

using RoME
using Test
using TensorCast
import Manifolds
using Manifolds: ProductManifold, SpecialEuclidean, ProductRepr, SpecialOrthogonal, TranslationGroup, identity_element
using DistributedFactorGraphs
using Statistics

# TODO move to IIF?
function Statistics.cov(vartype::InferenceVariable, ptsArr::AbstractVector; basis::Manifolds.AbstractBasis = Manifolds.DefaultOrthogonalBasis(), kwargs...)
  cov(getManifold(vartype), ptsArr; basis, kwargs... )
end


@error("must restore testG2oParser.jl")
@error("must restore testParametric.jl")
# "testG2oParser.jl";  ]

testfiles = [
  "testInflation380.jl";
  "testPoint2Point2.jl";
  # "testParametric.jl";
  "testTreeInitCommonMsg_IIF913.jl";
  "threeDimLinearProductTest.jl";
  "testPose3Pose3NH.jl";
  "testBeehive2D_CliqByCliq.jl";      # special case debugging
  "testhigherdimroots.jl";
  "testDidsonFunctions.jl";
  "testPoint2Point2Init.jl";
  "testBasicPose2Stationary.jl";
  "TestPoseAndPoint2Constraints.jl";
  "testPartialRangeCrossCorrelations.jl";
  "testDynPoint2D.jl";
  "testBearingRange2D.jl";
  "testBearing2D.jl";
  "testDeltaOdo.jl";
  "testFixedLagFG.jl";
  "testMultimodalRangeBearing.jl";
  "testDynPose2D.jl";
  "testPartialPriorYawPose2.jl";
  "testPartialPose3.jl";
  "testpackingconverters.jl";
  "TestDefaultFGInitialization.jl";
  "testAccumulateFactors.jl";
  "testDeadReckoningTether.jl"; 
  "testFluxModelsPose2.jl";
  "testBeehiveGrow.jl";
  "testGenerateHelix.jl";
]


## Tests not ready yet
# "HexagonalLightGraphs.jl"
# "testCameraFunctions.jl"
# "testmultiplefeatures.jl"


# test_results = @testset BrokenTestSet "Broken Testset for RoME" begin
for testf in testfiles
  println("[TEST] $testf")
  include(testf)
  println("[SUCCESS]")
  println()
  println()
  println()
end
# end

