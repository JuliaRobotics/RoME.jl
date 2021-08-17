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

## FIXME remove
#  quick and dirty skipping of broken testsets
mutable struct BrokenTestSet <: Test.AbstractTestSet
  description::AbstractString
  result::Vector
  n_pass::Int
  n_broken::Int
  n_fail::Int
  n_error::Int
end

BrokenTestSet(desc) = BrokenTestSet(desc,[],0,0,0,0)

Test.record(ts::BrokenTestSet, child::BrokenTestSet) = push!(ts.result, child)
Test.record(ts::BrokenTestSet, t::Test.Pass) = (ts.n_pass += 1; t)
Test.record(ts::BrokenTestSet, t::Test.Broken) = (ts.n_broken += 1; t)
Test.record(ts::BrokenTestSet, t::Test.Fail) = (ts.n_fail += 1; t)
Test.record(ts::BrokenTestSet, t::Test.Error) = (ts.n_error += 1; t)

function Test.finish(ts::BrokenTestSet)
  # just record if we're not the top-level parent
  if Test.get_testset_depth() > 0
      Test.record(Test.get_testset(), ts)
  end
  printstyled("Broken testset: ", ts.description, "\n"; bold=true)
  printstyled("Results:  ")
  printstyled("pass:$(ts.n_pass) ";color=:green)
  printstyled("broken:$(ts.n_broken) "; color=:yellow)
  printstyled("fail:$(ts.n_fail) "; color=:red)
  printstyled("error:$(ts.n_error)\n"; color=:red)
  ts
end
## end #FIXME remove

@error("must restore testG2oParser.jl")
# "testG2oParser.jl";  ]

testfiles = [
  "testInflation380.jl";
  "testPoint2Point2.jl";
  "testParametric.jl";
  "testTreeInitCommonMsg_IIF913.jl";
  "threeDimLinearProductTest.jl";
  "testPose3Pose3NH.jl";
  "testBeehive2D_CliqByCliq.jl";      # special case debugging
  "testhigherdimroots.jl";
  "testManifoldsPose2Equivalent.jl";
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
  "testPartialXYH.jl";
  "testpartialpose3.jl";
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
# "testPoint2Point2WorldBearing.jl";  # deprecate


test_results = @testset BrokenTestSet "Broken Testset for RoME" begin
for testf in testfiles
  println("[TEST] $testf")
  include(testf)
  println("[SUCCESS]")
  println()
  println()
  println()
end
end


n_pass = mapreduce(x->x.n_pass, +, test_results.result)
n_error = mapreduce(x->x.n_error, +, test_results.result)
n_fail = mapreduce(x->x.n_fail, +, test_results.result)
n_broken = mapreduce(x->x.n_broken, +, test_results.result)

printstyled("Broken Testset for RoME summary: \n"; bold=true)
printstyled("  Pass:$(n_pass)\n";color=:green)
printstyled("  Broken:$(n_broken)\n"; color=:yellow)
printstyled("  Fail:$(n_fail)\n"; color=:red)
printstyled("  Error:$(n_error)\n"; color=:red)

0 < n_error ? error("RoME tests failed") : nothing