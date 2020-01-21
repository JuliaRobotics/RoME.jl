# test Pose2Point2RangeNeighborhood factor

using Revise

using Test
using RoME
using IncrementalInference
using DistributedFactorGraphs



getMax(nn::Normal) = nn.μ

# must be the same as getSample
getMax(pprn::Union{Pose2Point2Range, Pose2Point2RangeNeighborhood}) = ([[getMax(pprn.Z)]';],[0.0;])

## Standard function which can act as parametric wrapper for factors
# kinda already have one

getMax(fct::FunctorPairwiseMinimize) = error("not implemented yet")
getMax(fct::FunctorPairwise) = error("not implemented yet")
getMax(fg::DFGFactor) = getMax(getFactorType(fc))


# return residual from the factor
function calcFactor(fcttype::FunctorPairwiseMinimize, vars...; measurement::Tuple=getMax(fc))
  zDim = size(measurement[1][1],2)
  res = zeros(zDim)
  ud = FactorMetadata()
  fcttype(res,ud,1,measurement,vars...)
  res
end

calcFactor(fc::DFGFactor, vars...; measurement::Tuple=getMax(fc)) = calFactor(getFactorType(fc), vars..., measurement=measurement)







@testset "Test Range2DNeighborhood factor residual function" begin


pprn = Pose2Point2RangeNeighborhood(Normal(10,0.3))


x1 = [0.0 0.0 0.0]'
l1 = [10.0 0.0]'

calcFactor(pprn, x1, l1, measurement=getMax(pprn))


x1 = [10.0 10.0 0.0]'

calcFactor(pprn, x1, l1, measurement=getMax(pprn))


x1 = [10.0 0.0 0.0]'

calcFactor(pprn, x1, l1, measurement=getMax(pprn))





end



@testset "testing Pose2Point2RangeNeighborhood" begin

fg = initfg()
getSolverParams(fg).N = 5

addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0;], PartialPriorYawPose2(Normal(0.0, 0.1) ), autoinit=false)

addVariable!(fg, :l0, Point2)
addVariable!(fg, :l1, Point2)
addFactor!(fg, [:l0], PriorPoint2(MvNormal([10.0; 0.0], [0.1 0; 0 0.1])) )
addFactor!(fg, [:l1], PriorPoint2(MvNormal([0.0; 10.0], [0.1 0; 0 0.1])) )


addFactor!(fg, [:x0;:l0], Pose2Point2RangeNeighborhood(Normal(10,0.3)) )
addFactor!(fg, [:x0;:l1], Pose2Point2RangeNeighborhood(Normal(10,0.3)) )

# drawGraph(fg, show=true)

ensureAllInitialized!(fg)


# existing values
getKDE(fg, :x0) |> getPoints


# project new values
# this result is not correct
pts = approxConv(fg, :x0l0f1, :x0)


manualinit!(fg, :x0, pts)


getKDE(fg, :l0) |> getPoints


# tree, smt, hist = solveTree!(fg)


using RoMEPlotting
Gadfly.set_default_plot_size(35cm, 25cm)


PP = manikde!(pts, Pose2)

plotPose(Pose2(), PP)




#
drawPosesLandms(fg)


reportFactors(fg, Pose2Point2RangeNeighborhood)


plotLocalProduct(fg, :x0)


plotPose(fg, :x0)









end




#
