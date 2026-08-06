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

# include("testDynPose2D.jl")

##

@error("add test for generateGraph_Beehive!, norm( simulated - default ) < tol")

testfiles = [
    # tests that run on AMP v0.15
    "testG2oParser.jl"
    "inertial/testIMUDeltaFactor.jl"
    "testPartialRangeCrossCorrelations.jl"
    "testParametric.jl"
    "testParametricCovariances.jl" # 2 test_broken
    "threeDimLinearProductTest.jl"
    "testBearingRange2D.jl"
    "testpackingconverters.jl" #FIXME for new DFG deprecations
    "testBasicPose2Conv.jl"
    "testhigherdimroots.jl"
    "testBasicPose2Conv.jl"
    "testDidsonFunctions.jl"
    "testDynPoint2D.jl"
    "testPose3.jl"
    "testDeltaOdo.jl"
    "TestDefaultFGInitialization.jl"
    "testAccumulateFactors.jl"
    "TestPoseAndPoint2Constraints.jl"
    "testBasicPose2Stationary.jl"
    "testPose2Propagate.jl"
    "testParametricSimulated.jl"
    "testGraphGenerators.jl"
    "testGenerateHelix.jl"


    "testPoint2Point2.jl" # numerics still a bit wild, which is less important during initial AMP v0.15 upgrade
    "testHexagonal2D_CliqByCliq.jl" 
    "testInflation380.jl"
    "testTreeInitCommonMsg_IIF913.jl"
    "testFixedLagFG.jl"
    "testDeadReckoningTether.jl" # SKIP accumulateMeans on solveParametric internal util step
    "testPoint2Point2Init.jl" # few skips


    "testPose3Pose3NH.jl" # slow
    "testDynPose2D.jl" # very slow, probably works, fixed L181, dim mismatch, expected len 5 got 8
]


#FIXME
fixme_broken_tests = [
    # known broken tests
    
    "inertial/testInertialDynamic.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    "inertial/testODE_INS.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    # "testFluxModelsPose2.jl";
    "testG2oExportSE3.jl" # FIX, The function `zero` exists, but no method is defined for this combination of argument types.
    
    "testVelPos3.jl" # FIX, NaN in json error
    # tests most likely to fail on numerics

    # recent development work
    "testPartialPose2.jl" # FIX, L39 numeric all ~zero
    "testPartialPose3.jl" # FIX, ManifoldPartials.jl:354 -- BoundsError: attempt to access 3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3) at index [3-element BitVector]
    "testBearing2D.jl" # FIX, BoundsError: attempt to access Float64 at index [2, 1], dim mismatch
    # "testMultimodalRangeBearing.jl" # restore after Bearing factors are fixed

    # regular tests expected to pass
    "testGenericProjection.jl" # broken COMPAT w CameraModels


    "testPartialPriorYawPose2.jl" # FIX partial [3] during HoDe bw optim

    "testScalarFields.jl" # SKIPPING SOLVE FOR SLOW AMP V0.15
    # broken PPE test
    "testBeehiveGrow.jl"
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
                include(testf)
        end
        println()
        println()
    end

    for (i, testf) in enumerate(fixme_broken_tests)
        @testset "[TEST $i] $testf" begin
            println(
                "[BROKEN TEST $i] $testf"
            )
            @test_broken false
        end
    end
end
