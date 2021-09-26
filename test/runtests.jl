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

testfiles = [
  # any wip fail-fast dev testing
  # "testPartialPose3.jl";
  
  # important tests are that are broken and must be restored.
  "testBearingRange2D.jl";
  "testBearing2D.jl";
  # "testMultimodalRangeBearing.jl"; # restore after Bearing factors are fixed
  
  # "testParametric.jl"; # deferred to v0.16.x
  # "testG2oParser.jl";  # deferred to v0.16.x
  
  # tests most likely to fail on numerics
  "testPoint2Point2Init.jl";
  "testBeehiveGrow.jl"; # also starts multiprocess
  
  # regular tests expected to pass
  "testpackingconverters.jl";
  "testInflation380.jl";
  "testPoint2Point2.jl";
  "testParametricSimulated.jl";
  "testBasicPose2Conv.jl";
  "testTreeInitCommonMsg_IIF913.jl";
  "threeDimLinearProductTest.jl";
  "testPose3Pose3NH.jl";
  "testHexagonal2D_CliqByCliq.jl";      # special case debugging
  "testhigherdimroots.jl";
  "testDidsonFunctions.jl";
  "testBasicPose2Stationary.jl";
  "TestPoseAndPoint2Constraints.jl";
  "testPartialRangeCrossCorrelations.jl";
  "testDynPoint2D.jl";
  "testDeltaOdo.jl";
  "testFixedLagFG.jl";
  "testAccumulateFactors.jl";
  "testDynPose2D.jl";
  "testPartialPriorYawPose2.jl";
  "TestDefaultFGInitialization.jl";
  "testDeadReckoningTether.jl"; 
  "testFluxModelsPose2.jl";
  "testGenerateHelix.jl";
]


## Tests not ready yet
# "HexagonalLightGraphs.jl"
# "testCameraFunctions.jl"
# "testmultiplefeatures.jl"


for testf in testfiles
  println("[TEST] $testf")
  include(testf)
  println("[SUCCESS] $testf ==========================================================")
  println()
  println()
  println()
end


#