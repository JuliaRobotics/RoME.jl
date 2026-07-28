#addprocs(2)
#@show nprocs()
using RoMETypes
using RoME
using Test
using TensorCast
import Manifolds
using LieGroups
using LieGroups: TranslationGroup
# using Manifolds: ProductManifold, SpecialEuclidean, SpecialOrthogonal, TranslationGroup, identity_element
using DistributedFactorGraphs
DFG.@usingDFG true
using Statistics
using LinearAlgebra
using Random
using StaticArrays

using RoME: LeftInvariantMetricSE
using ManifoldsBase: submanifold_component

##

# include("testScalarFields.jl")

##

@error("add test for generateGraph_Beehive!, norm( simulated - default ) < tol")

testfiles = [
    "inertial/testIMUDeltaFactor.jl"
    "testPartialRangeCrossCorrelations.jl"
    "testParametric.jl"
    
    # known broken tests
    "testG2oParser.jl"  # deferred
    "inertial/testInertialDynamic.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    "inertial/testODE_INS.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    # "testFluxModelsPose2.jl";
    "testG2oExportSE3.jl" # FIX
    "testPose2Propagate.jl" # FIX
    "testPose3.jl" # FIX
    "testVelPos3.jl" # FIX

    # tests most likely to fail on numerics
    "testScalarFields.jl" # SKIPPING SOLVE FOR SLOW AMP V0.15
    "testPoint2Point2Init.jl" # FIX
    "threeDimLinearProductTest.jl" # FIX
    "testPose3Pose3NH.jl" # FIX,   UndefVarError: `LeftInvariantMetricSE` not defined in `Main`

    # recent development work
    "testPartialPose2.jl" # FIX, L39 numeric all ~zero
    "testPartialPose3.jl" # FIX, UndefVarError: `LeftInvariantMetricSE` not defined in `Main` 
    "testBearingRange2D.jl"
    "testBearing2D.jl"
    "testMultimodalRangeBearing.jl" # restore after Bearing factors are fixed

    # regular tests expected to pass
    "testpackingconverters.jl" #FIXME for new DFG deprecations
    "testInflation380.jl"
    "testPoint2Point2.jl"
    "testParametricCovariances.jl"
    "testParametricSimulated.jl" # FIX
    "testBasicPose2Conv.jl"
    "testGraphGenerators.jl" # FIX
    "testTreeInitCommonMsg_IIF913.jl"
    "testHexagonal2D_CliqByCliq.jl"      # special case debugging
    "testhigherdimroots.jl"
    "testGenericProjection.jl"
    "testDidsonFunctions.jl"
    "testBasicPose2Stationary.jl"
    "TestPoseAndPoint2Constraints.jl"
    "testDynPoint2D.jl"
    "testDeltaOdo.jl"
    "testFixedLagFG.jl"
    "testDynPose2D.jl"
    "testPartialPriorYawPose2.jl"
    "TestDefaultFGInitialization.jl"
    "testAccumulateFactors.jl"
    "testDeadReckoningTether.jl"
    "testGenerateHelix.jl" # FIX

    # starts multiprocess.
    # don't move up, special factors defined in other test files are not added to multiprocess (Distributed.jl)
    # "testBeehiveGrow.jl" # FIX, also starts multiprocess
]

#FIXME
fixme_broken_ppe_tests = [
    "testScalarFields.jl",
    "testParametricSimulated.jl",
    "testGraphGenerators.jl",
    "testGenerateHelix.jl",
    "testBeehiveGrow.jl",
]

## Tests not ready yet
# "HexagonalLightGraphs.jl"
# "testCameraFunctions.jl"
# "testmultiplefeatures.jl"
@testset "RoME tests" begin
    for (i, testf) in enumerate(testfiles)
        @testset "[TEST $i] $testf" begin
            println(
                "[TEST $i] $testf =============================================================",
            )
            if testf in fixme_broken_ppe_tests
                @test_broken false
            else
                include(testf)
            end
        end
        println()
        println()
        println()
    end
end
