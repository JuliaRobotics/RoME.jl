# test Pose2DPoint2D constraint evaluation function

using RoME, IncrementalInference
using Base.Test

# begin
N = 100
fg = initfg()


initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.01])

# Some starting position
v1 = addNode!(fg, :x1, zeros(3,1), diagm([1.0;1.0;0.1]), N=N)
# and a second pose, note there is now prior to give gauge
v2 = addNode!(fg, :x2, ([50.0;0.0;0.0]')', diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;0.0]')', 0.1*eye(3), [1.0])
f2 = addFactor!(fg, [v1;v2], ppc)

# test evaluation of pose pose constraint
pts = evalFactor2(fg, f2, v2.index)
@test norm(Base.mean(pts,2)[1:2]-[50.0;0.0]) < 10.0
@test abs(Base.mean(pts,2)[3]) < 0.5

# new landmark
l1 = addNode!(fg, :l1, ([10.0;0.0]')', diagm([1.0;1.0]), N=N)
# and pose to landmark constraint
rhoZ1 = norm([10.0;0.0])
ppr = Pose2DPoint2DRange([rhoZ1], 2.0, [1.0])
f1 = addFactor!(fg, [v1;l1], ppr)

pts = evalFactor2(fg, f1, l1.index)

# using KernelDensityEstimate
# using Gadfly
#
# p1= kde!(pts)
#
# plotKDE( marginal(p1,[1;2]) , dimLbls=["x";"y";"z"])
# plotKDE( marginal(getVertKDE(fg, :x2),[1;2]) , dimLbls=["x";"y";"z"])
#
# axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )
#

# @show sum(sqrt(sum(pts.^2, 1 )) .< 10.0)
# @test sum(sqrt(sum(pts.^2, 1 )) .< 10.0) == 0



# add a prior somewhere
pp2 = PriorPoint2D([10.0;0.0], diagm([1.0;1.0]), [1.0])

f3 = addFactor!(fg,[l1], pp2)
pts = evalFactor2(fg, f3, l1.index)

@test norm(Base.mean(pts,2)[:]-[10.0;0.0]) < 3.0
# end
