#addprocs(2)
#@show nprocs()

using RoME
using Test
using Pkg


pkgversion(m::Module) = Pkg.TOML.parsefile(joinpath(dirname(string(first(methods(m.eval)).file)), "..", "Project.toml"))["version"]

@show pkgversion(IncrementalInference)

function _defaultFactorMetadataRoME(Xi::AbstractVector{<:DFGVariable};
                                solvefor::Symbol=:null,
                                arrRef=Vector{Matrix{Float64}}(),
                                cachedata::T=nothing ) where T
  #
  

  #NOTE temp for IIF v0.20 migration
  # FIXME standardize fmd, see #927
  iifv = pkgversion(IncrementalInference) |> string |> VersionNumber
  ret = if iifv < v"0.20"
    FactorMetadata(solvefor,map(x->x.label,Xi),cachedata,copy(Xi), arrRef)
  else
    FactorMetadata(copy(Xi), map(x->x.label,Xi), arrRef, solvefor, cachedata)
  end
  return ret
end

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
"testDeadReckoningTether.jl"; ]

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
