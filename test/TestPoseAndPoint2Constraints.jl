# test Pose2DPoint2D constraint evaluation function

using LinearAlgebra, Statistics
using Test
using RoME
# , IncrementalInference, Distributions

##

@testset "test pose and point combinations..." begin

##

N = 100
fg = initfg()


initCov = Matrix(Diagonal([0.03;0.03;0.001]))
odoCov = Matrix(Diagonal([3.0;3.0;0.01]))

# Some starting position
v1 = addVariable!(fg, :x0, Pose2, N=N)
initPosePrior = PriorPose2(MvNormal(zeros(3), initCov))
f1  = addFactor!(fg,[v1], initPosePrior)

@test Pose2Pose2(MvNormal(randn(2), Matrix{Float64}(LinearAlgebra.I, 2,2))) !== nothing
@test Pose2Pose2(MvNormal(randn(2), Matrix{Float64}(LinearAlgebra.I, 2,2))) !== nothing
@test Pose2Pose2(MvNormal(randn(2), Matrix{Float64}(LinearAlgebra.I, 2,2))) !== nothing

# and a second pose
v2 = addVariable!(fg, :x1, Pose2, N=N)
ppc = Pose2Pose2(MvNormal([50.0;0.0;pi/2], odoCov))
f2 = addFactor!(fg, [:x0;:x1], ppc)

# test evaluation of pose pose constraint
pts = approxConv(fg, :x0x1f1, :x1)

# pts[3,:] .= TU.wrapRad.(pts[3,:])
M = getManifold(Pose2)
@show mv = mean(M, pts)
@test isapprox(submanifold_component(mv,1), [50,0], atol=1)
@test isapprox(submanifold_component(mv,2), [0 -1; 1 0], atol=0.5)

##

# @show ls(fg)
tree = solveTree!(fg)
# initAll!(fg)
# tree = wipeBuildNewTree!(fg)
# inferOverTreeR!(fg, tree,N=N)
# inferOverTree!(fg, tree, N=N)

# test post evaluation values are correct
pts = getVal(fg, :x0)
@test isapprox(M, mean(M, pts), ArrayPartition([0,0], [1 0; 0 1]), atol=0.5)


pts = getVal(fg, :x1)
me_ = mean(M, pts)
@test isapprox(M.manifold[1], submanifold_component(me_,1), [50,0], atol=1.0)
@test isapprox(M.manifold[2], submanifold_component(me_,2), [0 -1; 1 0], atol=0.25)

# check that yaw is working
v3 = addVariable!(fg, :x2, Pose2, N=N)
ppc = Pose2Pose2(MvNormal([50.0;0.0;0.0], odoCov))
f3 = addFactor!(fg, [v2;v3], ppc)

solveTree!(fg)
# solveTree!(fg)

# test post evaluation values are correct
pts = getVal(fg, :x0)
me_ = mean(M, pts)
@test isapprox(M.manifold[1], submanifold_component(me_,1), [0,0], atol=0.5)
@test isapprox(M.manifold[2], submanifold_component(me_,2), [1 0; 0 1], atol=0.25)

pts = getVal(fg, :x1)
mv = mean(M, pts)
@test isapprox(submanifold_component(mv,1), [50,0], atol=3.0)
@test isapprox(submanifold_component(mv,2), [0 -1; 1 0], atol=0.1)

pts = getVal(fg, :x2)
mv = mean(M, pts)
@test isapprox(submanifold_component(mv,1), [50,50], atol=3.0)
@test isapprox(submanifold_component(mv,2), [0 -1; 1 0], atol=0.1)

println("test bearing range evaluations")

# new landmark
l1 = addVariable!(fg, :l1, Point2, N=N)
# and pose to landmark constraint
rhoZ1 = norm([10.0;0.0])
ppr = Pose2Point2BearingRange(Uniform(-pi,pi),Normal(rhoZ1,1.0))
f4 = addFactor!(fg, [:x0;:l1], ppr)


pts = approxConv(fg, :x0l1f1, :l1, N=N)
# pts = evalFactor(fg, f4, l1.index, N=N)
## res = zeros(2)
## meas = getSample(ppr)
## ppr(res,nothing,1,meas,xi,lm )

# all points should lie in a ring around 0,0
@test sum(norm.(pts) .< 5.0) == 0
@test sum(norm.(pts) .< 15.0) == N


pts = approxConv(fg, :x0l1f1, :x0)
# pts = evalFactor(fg, f4, v1.index, N=200)
# @show sum(sqrt(sum(pts.^2, 1 )) .< 5.0)
# @test sum(sqrt(sum(pts.^2, 1 )) .< 5.0) == 0



# add a prior to landmark
pp2 = PriorPoint2(MvNormal([10.0;0.0], Matrix(Diagonal([1.0;1.0]))))

f5 = addFactor!(fg,[l1], pp2)
pts = approxConv(fg, f5.label, l1.label)
# pts = evalFactor(fg, f5, l1.label)


@test isapprox(Statistics.mean(pts), [10, 0], atol=1.0)

# println("test Pose2D plotting")

# drawPoses(fg);
# drawPosesLandms(fg);

# using KernelDensityEstimate
# using Gadfly
#
# pts = getVal(fg, :l1)
#
# p1= kde!(pts)
# plotKDE( p1 , dimLbls=["x";"y";"z"])

# plotKDE( [marginal(p1c,[1;2]);marginal(p1,[1;2])] , dimLbls=["x";"y";"z"],c=["red";"black"],levels=3)
# p1c = deepcopy(p1)

# plotKDE( marginal(getBelief(fg, :x2),[1;2]) , dimLbls=["x";"y";"z"])
#
# axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )
#


end




#
