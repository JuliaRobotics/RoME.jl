#  test MickeyMouse2D
using RoME, Distributions
using IncrementalInference   # for evalFactor2
using Base.Test

# using Gadfly

x1 = [0.0;0;0]
x2 = [1.0;-1.0;pi/2]

l1 = [0.5;0.5]
l2 = [1.0;0.0]
l3 = [1.5;-0.5]
L = [l1';l2';l3']

# plot(x=L[:,1],y=L[:,2],Geom.point)

the11 = atan2(l1[2],l1[1])
the12 = atan2(l2[2],l2[1])
the13 = atan2(l3[2],l3[1])

Xi = SE2(x1)
Xj = SE2(x2)
Δxij = se2vee(Xi\Xj)

tl1 = se2vee(Xj \ SE2([l1;0.0]))[1:2]
the21 = atan2(tl1[2],tl1[1])

tl2 = se2vee(Xj \ SE2([l2;0.0]))[1:2]
the22 = atan2(tl2[2],tl2[1])

tl3 = se2vee(Xj \ SE2([l3;0.0]))[1:2]
the23 = atan2(tl3[2],tl3[1])


xir1 = Normal(the11,1e-4)
xir2 = Normal(the12,1e-4)
xir3 = Normal(the13,1e-4)
xjr1 = Normal(the21,1e-4)
xjr2 = Normal(the22,1e-4)
xjr3 = Normal(the23,1e-4)

mm2 = MickeyMouse2D(
  xir1,
  xir2,
  xir3,
  xjr1,
  xjr2,
  xjr3,
  SE2(zeros(3))  )

res = randn(2)
idx = 1
meas = getSample(mm2)


wAbi = (x1')'
wAbj = (x2')'
wAo1 = ([l1;0.0]')'
wAo2 = ([l2;0.0]')'
wAo3 = ([l3;0.0]')'


mm2(res, idx, meas,
  wAbi,
  wAbj,
  wAo1,
  wAo2,
  wAo3  )


@test norm(res) < 1e-3



# now build in factor graph form for further testing

N=100
fg = initfg()


v1 = addNode!(fg, :x1, 0.001*randn(3,N), diagm([1.0;1.0;0.1]), N=N)

pts2 = [0.1*randn(1,N)+1;  0.1*randn(1,N)-1; 0.01*randn(1,N)+pi/2]
v2 = addNode!(fg, :x2, pts2, diagm([1.0;1.0;0.1]), N=N)


vl1 = addNode!(fg, :l1, randn(2,N), diagm([1.0;1.0]), N=N)
vl2 = addNode!(fg, :l2, randn(2,N), diagm([1.0;1.0]), N=N)
vl3 = addNode!(fg, :l3, randn(2,N), diagm([1.0;1.0]), N=N)


f1  = addFactor!(fg, [v1;v2;vl1;vl2;vl3], mm2)


priorpts = evalFactor2(fg, f1, fg.IDs[:l1])




























#
