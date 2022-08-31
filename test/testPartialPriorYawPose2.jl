using RoME
using Test
using Manifolds: ProductRepr

##

@testset "Basic PartialPriorYawPose2 test" begin
##

fg = initfg()

# fg.solverParams.graphinit=false
addVariable!(fg, :x, Pose2)
addVariable!(fg, :l, Point2)

addFactor!(fg, [:x], PartialPriorYawPose2(MvNormal([0.,0, 0], [0.01, 0.01, 0.01])))
addFactor!(fg, [:l], PriorPoint2(MvNormal([1.,1], [0.01, 0.01])))

addFactor!(fg, [:x; :l], Pose2Point2(MvNormal([1.,1], [0.1,0.1])))

initAll!(fg)

tree = solveTree!(fg)

M = getManifold(Pose2)
me_ = mean(M, getVal(fg, :x))
#x should form a dounut around 1,1 with yaw close to zero
@test isapprox(M.manifold[1], submanifold_component(me_,1), [1,1], atol=0.2)
@test isapprox(M.manifold[2], submanifold_component(me_,2), [1 0; 0 1], atol=0.05) 
@test isapprox(mean(getVal(fg, :l)), [1,1], atol = 0.05)

##
end


@testset "Basic PartialPriorYawPose2 test with Bearing Range" begin
##

fg = initfg()

# fg.solverParams.graphinit=false
addVariable!(fg, :x, Pose2)
addVariable!(fg, :l, Point2)

addFactor!(fg, [:x], PartialPriorYawPose2(Normal(0.0, 0.001)))
addFactor!(fg, [:l], PriorPoint2(MvNormal([10.,10], [0.01, 0.01])))

addFactor!(fg, [:x; :l], Pose2Point2BearingRange(Normal(0.0, 0.01),Normal(5.0, 0.1)))

initAll!(fg)

tree = solveTree!(fg)

M = getManifold(Pose2)
me_ = mean(M, getVal(fg, :x))
@test isapprox(M.manifold[1], submanifold_component(me_,1), [5,10], atol=0.2)
@test isapprox(M.manifold[2], submanifold_component(me_,2), [1 0; 0 1], atol=0.05) 
@test isapprox(mean(getVal(fg, :l)), [10,10], atol = 0.05)

# ##
end

# Debug plots using Makie
# points = getPoints(fg, :x)
# @cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
# f = scatter(pts[1,:], pts[2,:]; rotations = pts[3,:], markersize=15, marker = '►')
# xlims!(0,20)
# ylims!(0,20)
# f.axis.aspect[] =1.0
# points = getPoints(fg, :l)
# @cast pts[j,i] := getCoordinates.(Point2,points)[i][j]
# scatter!(pts[1,:], pts[2,:])

##