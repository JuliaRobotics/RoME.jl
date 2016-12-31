# test Pose2DPoint2D constraint evaluation function

using RoME, IncrementalInference, Distributions
using Base.Test

begin
N = 100
fg = initfg()


initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.01])

# Some starting position
v1 = addNode!(fg, :x1, zeros(3,1), diagm([1.0;1.0;0.1]), N=N)
initPosePrior = PriorPose2(zeros(3,1), initCov, [1.0])
f1  = addFactor!(fg,[v1], initPosePrior)

# and a second pose
v2 = addNode!(fg, :x2, ([50.0;0.0;pi/2]')', diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;pi/2]')', odoCov, [1.0])
f2 = addFactor!(fg, [v1;v2], ppc)

# test evaluation of pose pose constraint
pts = evalFactor2(fg, f2, v2.index)
@test norm(Base.mean(pts,2)[1:2]-[50.0;0.0]) < 10.0
@test abs(Base.mean(pts,2)[3]-pi/2) < 0.5

# @show ls(fg)

tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree)
inferOverTree!(fg, tree)

# test post evaluation values are correct
pts = getVal(fg, :x1)
@test norm(Base.mean(pts,2)[1:2]-[0.0;0.0]) < 10.0
@test abs(Base.mean(pts,2)[3]) < 0.5

pts = getVal(fg, :x2)
@test norm(Base.mean(pts,2)[1:2]-[50.0;0.0]) < 10.0
@test abs(Base.mean(pts,2)[3]-pi/2) < 0.5


# check that yaw is working
v3 = addNode!(fg, :x3, ([0.0;0.0;0.0]')', diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;0.0]')', odoCov, [1.0])
f3 = addFactor!(fg, [v2;v3], ppc)


tree = wipeBuildNewTree!(fg)
[inferOverTree!(fg, tree) for i in 1:3]

# test post evaluation values are correct
pts = getVal(fg, :x1)
@test norm(Base.mean(pts,2)[1:2]-[0.0;0.0]) < 20.0
@test abs(Base.mean(pts,2)[3]) < 0.5

pts = getVal(fg, :x2)
@test norm(Base.mean(pts,2)[1:2]-[50.0;0.0]) < 20.0
@test abs(Base.mean(pts,2)[3] - pi/2) < 0.5

pts = getVal(fg, :x3)
@test norm(Base.mean(pts,2)[1:2]-[50.0;50.0]) < 20.0
@test abs(Base.mean(pts,2)[3]-pi/2) < 0.5


# new landmark
l1 = addNode!(fg, :l1, ([0.0;0.0]')', diagm([1.0;1.0]), N=N)
# and pose to landmark constraint
rhoZ1 = norm([10.0;0.0])
ppr = Pose2DPoint2DBearingRange{Uniform, Normal}(Uniform(-pi,pi),Normal(rhoZ1,1.0))
f4 = addFactor!(fg, [v1;l1], ppr)

pts = evalFactor2(fg, f4, l1.index)
# @show sum(sqrt(sum(pts.^2, 1 )) .< 5.0)
@test sum(sqrt(sum(pts.^2, 1 )) .< 5.0) == 0



# add a prior somewhere
pp2 = PriorPoint2D([10.0;0.0], diagm([1.0;1.0]), [1.0])

f5 = addFactor!(fg,[l1], pp2)
pts = evalFactor2(fg, f5, l1.index)

@test norm(Base.mean(pts,2)[:]-[10.0;0.0]) < 5.0


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
