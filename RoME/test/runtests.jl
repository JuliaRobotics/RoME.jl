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

# include("testPose3.jl")

##

@error("add test for generateGraph_Beehive!, norm( simulated - default ) < tol")

testfiles = [
    # tests that run on AMP v0.15
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

    "testPoint2Point2.jl" # numerics still a bit wild, which is less important during initial AMP v0.15 upgrade


    # known broken tests
    "testG2oParser.jl"  # deferred
    "inertial/testInertialDynamic.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    "inertial/testODE_INS.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    # "testFluxModelsPose2.jl";
    "testG2oExportSE3.jl" # FIX
    "testPose2Propagate.jl" # FIX
    
    "testVelPos3.jl" # FIX

    # tests most likely to fail on numerics
    "testScalarFields.jl" # SKIPPING SOLVE FOR SLOW AMP V0.15
    "testPoint2Point2Init.jl" # FIX
    "testPose3Pose3NH.jl" # FIX, dim mismatch, expected len 6 got 12

    # recent development work
    "testPartialPose2.jl" # FIX, L39 numeric all ~zero
    "testPartialPose3.jl" # FIX, MethodError: no method matching getManifoldPartial(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, ::Vector{Int64}, ::SMatrix{3, 3, Float64, 9}, ::Base.RefValue{Int64}; doError::Bool)
    "testBearing2D.jl" # FIX, dim mismatch expected len 3 got 6
    # "testMultimodalRangeBearing.jl" # restore after Bearing factors are fixed

    # regular tests expected to pass
    "testInflation380.jl" # FIX, expected input len 3 got 6

    "testParametricSimulated.jl" # FIX
    "testGraphGenerators.jl" # FIX
    "testTreeInitCommonMsg_IIF913.jl" # FIX dim mismatch, expected len 3 got 6
    "testHexagonal2D_CliqByCliq.jl"  # FIX, dim mismatch, expected len 3 got 6 # special case debugging

    "testGenericProjection.jl" # broken COMPAT w CameraModels
    "testBasicPose2Stationary.jl" # FIX

    "TestPoseAndPoint2Constraints.jl" # FIX, dim mismatch, expected len 3 got 6

    "testDynPose2D.jl" # FIX, L181, dim mismatch, expected len 5 got 8
    "testFixedLagFG.jl" # FIX, dim mismatch, expected len 3 got 6
    "testPartialPriorYawPose2.jl" # FIX partial [3] during HoDe bw optim
    "testDeadReckoningTether.jl" # FIX, dim mismatch, expected len 3 got 6
    "testGenerateHelix.jl" # FIX

    # starts multiprocess.
    # don't move up, special factors defined in other test files are not added to multiprocess (Distributed.jl)
    # "testBeehiveGrow.jl" # FIX, also starts multiprocess
]

#FIXME
fixme_broken_tests = [
    # known broken tests
    "testG2oParser.jl"  # deferred
    "inertial/testInertialDynamic.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    "inertial/testODE_INS.jl" # same problem as IIF/test/testDERelative.jl, access [1] of [0]
    # "testFluxModelsPose2.jl";
    "testG2oExportSE3.jl" # FIX
    "testPose2Propagate.jl" # FIX
    "testVelPos3.jl" # FIX
    # tests most likely to fail on numerics
    "testScalarFields.jl" # SKIPPING SOLVE FOR SLOW AMP V0.15
    "testPoint2Point2Init.jl" # FIX
    "testPose3Pose3NH.jl" # FIX, dim mismatch, expected len 6 got 12
    # recent development work
    "testPartialPose2.jl" # FIX, L39 numeric all ~zero
    "testPartialPose3.jl" # FIX, MethodError: no method matching getManifoldPartial(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, ::Vector{Int64}, ::SMatrix{3, 3, Float64, 9}, ::Base.RefValue{Int64}; doError::Bool)
    "testBearing2D.jl" # FIX, dim mismatch expected len 3 got 6
    # "testMultimodalRangeBearing.jl" # restore after Bearing factors are fixed
    # regular tests expected to pass
    "testInflation380.jl" # FIX, expected input len 3 got 6
    "testParametricSimulated.jl" # FIX
    "testGraphGenerators.jl" # FIX
    "testTreeInitCommonMsg_IIF913.jl" # FIX dim mismatch, expected len 3 got 6
    "testHexagonal2D_CliqByCliq.jl"  # FIX, dim mismatch, expected len 3 got 6 # special case debugging
    "testGenericProjection.jl" # broken COMPAT w CameraModels
    "testBasicPose2Stationary.jl" # FIX
    "TestPoseAndPoint2Constraints.jl" # FIX, dim mismatch, expected len 3 got 6
    "testDynPose2D.jl" # FIX, L181, dim mismatch, expected len 5 got 8
    "testFixedLagFG.jl" # FIX, dim mismatch, expected len 3 got 6
    "testPartialPriorYawPose2.jl" # FIX partial [3] during HoDe bw optim
    "testDeadReckoningTether.jl" # FIX, dim mismatch, expected len 3 got 6
    "testGenerateHelix.jl" # FIX


    # broken PPE tests
    "testScalarFields.jl"
    "testParametricSimulated.jl"
    "testGraphGenerators.jl"
    "testGenerateHelix.jl"
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
            if testf in fixme_broken_tests
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
