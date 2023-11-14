#addprocs(2)
#@show nprocs()

using RoME
using Test
using TensorCast
import Manifolds
using Manifolds: ProductManifold, SpecialEuclidean, SpecialOrthogonal, TranslationGroup, identity_element
using DistributedFactorGraphs
using Statistics
using StaticArrays


@error("add test for generateGraph_Beehive!, norm( simulated - default ) < tol")

testfiles = [  
  # known broken tests
  "testG2oParser.jl";  # deferred
  
  # dev test first, for faster issues.
  # Inertial
  "inertial/testODE_INS.jl";
  "inertial/testIMUDeltaFactor.jl";
  "inertial/testInertialDynamic.jl";
  
  # ...
  # "testFluxModelsPose2.jl";
  "testPartialRangeCrossCorrelations.jl";
  "testG2oExportSE3.jl";

  #parametric tests
  "testParametric.jl";
  "testPose3.jl";

  # tests most likely to fail on numerics
  "testScalarFields.jl";
  "testPoint2Point2Init.jl";
  "threeDimLinearProductTest.jl";
  "testPose3Pose3NH.jl";

  # recent development work
  "testPartialPose2.jl";
  "testPartialPose3.jl";
  "testBearingRange2D.jl";
  "testBearing2D.jl";
  "testMultimodalRangeBearing.jl"; # restore after Bearing factors are fixed

  # regular tests expected to pass
  "testpackingconverters.jl";
  "testInflation380.jl";
  "testPoint2Point2.jl";
  "testParametricCovariances.jl";
  "testParametricSimulated.jl";
  "testBasicPose2Conv.jl";
  "testGraphGenerators.jl";
  "testTreeInitCommonMsg_IIF913.jl";
  "testHexagonal2D_CliqByCliq.jl";      # special case debugging
  "testhigherdimroots.jl";
  "testGenericProjection.jl";
  "testDidsonFunctions.jl";
  "testBasicPose2Stationary.jl";
  "TestPoseAndPoint2Constraints.jl";
  "testDynPoint2D.jl";
  "testDeltaOdo.jl";
  "testFixedLagFG.jl";
  "testDynPose2D.jl";
  "testPartialPriorYawPose2.jl";
  "TestDefaultFGInitialization.jl";
  "testAccumulateFactors.jl";
  "testDeadReckoningTether.jl"; 
  "testGenerateHelix.jl";


  # starts multiprocess.
  # don't move up, special factors defined in other test files are not added to multiprocess (Distributed.jl)
  "testBeehiveGrow.jl"; # also starts multiprocess
]


## Tests not ready yet
# "HexagonalLightGraphs.jl"
# "testCameraFunctions.jl"
# "testmultiplefeatures.jl"

for (i,testf) in enumerate(testfiles)
  println("[TEST $i] $testf =============================================================")
  include(testf)
  println("[SUCCESS] $testf")
  println()
  println()
  println()
end


#