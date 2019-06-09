# test Pose2DPoint2D constraint evaluation function

using LinearAlgebra, Statistics
using Test
using RoME
# , IncrementalInference, Distributions


@testset "test pose and point combinations..." begin


N = 100
fg = initfg()


initCov = Matrix(Diagonal([0.03;0.03;0.001]))
odoCov = Matrix(Diagonal([3.0;3.0;0.01]))

# Some starting position
v1 = addVariable!(fg, :x0, Pose2, N=N)
initPosePrior = PriorPose2(MvNormal(zeros(3), initCov))
f1  = addFactor!(fg,[v1], initPosePrior)

@test Pose2Pose2(MvNormal(randn(2), Matrix{Float64}(LinearAlgebra.I, 2,2))) != nothing
@test Pose2Pose2(randn(2), Matrix{Float64}(LinearAlgebra.I, 2,2)) != nothing
@test Pose2Pose2(randn(2), Matrix{Float64}(LinearAlgebra.I, 2,2),[1.0;]) != nothing

# and a second pose
v2 = addVariable!(fg, :x1, Pose2, N=N)
ppc = Pose2Pose2([50.0;0.0;pi/2], odoCov, [1.0])
f2 = addFactor!(fg, [:x0;:x1], ppc)

# test evaluation of pose pose constraint
pts = approxConv(fg, :x0x1f1, :x1)
# pts = evalFactor2(fg, f2, v2.index)
@test norm(Statistics.mean(pts,dims=2)[1:2] - [50.0;0.0]) < 10.0
@test abs(Statistics.mean(pts,dims=2)[3] - pi/2) < 0.5

# @show ls(fg)
tree = batchSolve!(fg, treeinit=true)
# ensureAllInitialized!(fg)
# tree = wipeBuildNewTree!(fg)
# inferOverTreeR!(fg, tree,N=N)
# inferOverTree!(fg, tree, N=N)

# test post evaluation values are correct
pts = getVal(fg, :x0)
@test norm(Statistics.mean(pts, dims=2)[1:2] - [0.0;0.0]) < 10.0
@test abs(Statistics.mean(pts, dims=2)[3]) < 0.5

pts = getVal(fg, :x1)
@test norm(Statistics.mean(pts, dims=2)[1:2]-[50.0;0.0]) < 10.0
@test abs(Statistics.mean(pts, dims=2)[3]-pi/2) < 0.5


# check that yaw is working
v3 = addVariable!(fg, :x2, Pose2, N=N)
ppc = Pose2Pose2([50.0;0.0;0.0], odoCov, [1.0])
f3 = addFactor!(fg, [v2;v3], ppc)

batchSolve!(fg, treeinit=true)
batchSolve!(fg, treeinit=true)
# tree = wipeBuildNewTree!(fg)
# [inferOverTree!(fg, tree, N=N) for i in 1:2]

# test post evaluation values are correct
pts = getVal(fg, :x0)
@test norm(Statistics.mean(pts, dims=2)[1:2]-[0.0;0.0]) < 20.0
@test abs(Statistics.mean(pts, dims=2)[3]) < 0.5

pts = getVal(fg, :x1)
@test norm(Statistics.mean(pts, dims=2)[1:2]-[50.0;0.0]) < 20.0
@test abs(Statistics.mean(pts, dims=2)[3] - pi/2) < 0.5

pts = getVal(fg, :x2)
@test norm(Statistics.mean(pts, dims=2)[1:2]-[50.0;50.0]) < 20.0
@test abs(Statistics.mean(pts, dims=2)[3]-pi/2) < 0.5

println("test bearing range evaluations")

# new landmark
l1 = addVariable!(fg, :l1, Point2, N=N)
# and pose to landmark constraint
rhoZ1 = norm([10.0;0.0])
ppr = Pose2Point2BearingRange(Uniform(-pi,pi),Normal(rhoZ1,1.0))
f4 = addFactor!(fg, [:x0;:l1], ppr)


pts = approxConv(fg, :x0l1f1, :l1, N=N)
# pts = evalFactor2(fg, f4, l1.index, N=N)
## res = zeros(2)
## meas = getSample(ppr)
## ppr(res,nothing,1,meas,xi,lm )

# all points should lie in a ring around 0,0
@test sum(sqrt.(sum(pts.^2, dims=1 )) .< 5.0) == 0
@test sum(sqrt.(sum(pts.^2, dims=1 )) .< 15.0) == N


pts = approxConv(fg, :x0l1f1, :x0)
# pts = evalFactor2(fg, f4, v1.index, N=200)
# @show sum(sqrt(sum(pts.^2, 1 )) .< 5.0)
# @test sum(sqrt(sum(pts.^2, 1 )) .< 5.0) == 0



# add a prior to landmark
pp2 = PriorPoint2(MvNormal([10.0;0.0], Matrix(Diagonal([1.0;1.0]))))

f5 = addFactor!(fg,[l1], pp2)
pts = approxConv(fg, f5.label, l1.label)
# pts = evalFactor2(fg, f5, l1.label)


@test norm(Statistics.mean(pts, dims=2)[:]-[10.0;0.0]) < 5.0

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

# plotKDE( marginal(getVertKDE(fg, :x2),[1;2]) , dimLbls=["x";"y";"z"])
#
# axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )
#


end




#
