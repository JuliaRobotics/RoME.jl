#  test MickeyMouse2D
using RoME, Distributions
using IncrementalInference   # for evalFactor2
using Base.Test

using Optim

using Gadfly

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

res = randn(3)
idx = 1
meas = getSample(mm2)


wAbi = (x1')'
wAbj = (x2')'
wAo1 = (l1')'
wAo2 = (l2')'
wAo3 = (l3')'


println("test MickeyMouse unit vector functions")
res, rhat, resid = zeros(3),zeros(2), zeros(2)

for i in 1:3
  r1, rhathat1, α1 = getUvecScaleFeature2D(Xi, SE2(zeros(3)), SE2([L[i,:];0.0]))
  rhat[1] = cos(meas[1][i,idx,1])
  rhat[2] = sin(meas[1][i,idx,1])
  @show rhat, rhathat1
  resid[1:2] = rhat-rhathat1
  @show res[1] += dot(resid, resid)
end

@test 0.0 <= res[1] < 1e-5


println("test residual function")

res = zeros(3)

wAbi[:] = 0.0
# @show wAbi

mm2(res, idx, meas,
  wAbi,
  wAbj,
  wAo1,
  wAo2,
  wAo3  )

# @show res

function gg2(x::Float64,y::Float64)
  wAbi[2] = x
  wAbi[3] = y

  res = zeros(3)

  mm2(res, idx, meas,
    wAbi,
    wAbj,
    wAo1,
    wAo2,
    wAo3  )
  return sum(res)[1]
end
XX = linspace(-0.5,0.5,200)
YY = linspace(-pi/4,pi/4,100)
gg2(0.0,0.0)
# Gadfly.plot(z=gg2,x=XX,y=YY, Geom.contour)

# ar99
println("test as minimization problem to pose")
function minmickey(x::Array{Float64,1})
  wAbi[1:3,1] = x[1:3]
  res = zeros(3)
  mm2(res, idx, meas,
    wAbi,
    wAbj,
    wAo1,
    wAo2,
    wAo3  )

  return sum(res[1:2])
end


gg = (x,y) -> minmickey([x;y;0.0])
# gg = (x,y) -> x*y
# plot(z=gg, x=linspace(-2,2,100),y=linspace(-2,2,100), Geom.contour)


gg(0.0,0.1)

@time r = optimize(minmickey, randn(3))
@show r.minimizer

@test norm(res) < 1e-3


println("test as minimization problem to landmark")
wAbi[:] = 0.0
function minmickey(x::Array{Float64,1})
  wAo1[1:2,1] = x[1:2]
  res = zeros(3)
  mm2(res, idx, meas,
    wAbi,
    wAbj,
    wAo1,
    wAo2,
    wAo3  )

  return sum(res[1:2])
end

gg = (x,y) -> minmickey([x;y])
#plot(z=gg, x=linspace(1,2,100),y=linspace(-1,1,100), Geom.contour(levels=50))

# now build in factor graph form for further testing

N=100
fg = initfg()

initCov = 0.1*eye(3)
initCov[3,3] = 0.01
odoCov = deepcopy(initCov)

ipp = PriorPose2(zeros(3,1), initCov,[1.0])

v1 = addNode!(fg, :x1, 0.001*randn(3,N), diagm([1.0;1.0;0.1]), N=N)

f1  = addFactor!(fg,[:x1], ipp)


pts2 = [0.1*randn(1,N)+1;  0.1*randn(1,N)-1; 0.01*randn(1,N)+pi/2]
v2 = addNode!(fg, :x2, pts2, diagm([1.0;1.0;0.1]), N=N)

px2 = zeros(3,1)
px2[1:3,1] = [1.0;-1.0;pi/2]
ipp2 = PriorPose2(px2, initCov,[1.0])
f1  = addFactor!(fg,[:x2], ipp2)


vl1 = addNode!(fg, :l1, rand(MvNormal(l1,0.001*eye(2)),N), diagm([1.0;1.0]), N=N)
vl2 = addNode!(fg, :l2, rand(MvNormal(l2,0.001*eye(2)),N), diagm([1.0;1.0]), N=N)
vl3 = addNode!(fg, :l3, rand(MvNormal(l3,0.001*eye(2)),N), diagm([1.0;1.0]), N=N)


f2 = addFactor!(fg, [v1;v2;vl1;vl2;vl3], mm2)

ef2pts = evalFactor2(fg, f2, fg.IDs[:x2])

# plotKDE(kde!(ef2pts[1:2,:]))
# plotKDE(marginal(kde!(ef2pts),[1;2]))

tree = wipeBuildNewTree!(fg)
# spyCliqMat(tree.cliques[1])

inferOverTreeR!(fg,tree)
#
# using KernelDensityEstimate
#
#
# plotKDE([getVertKDE(fg, :l1);getVertKDE(fg, :l2);getVertKDE(fg, :l3)])
#
# plotKDE([marginal(getVertKDE(fg, :x1),[1;2]);marginal(getVertKDE(fg, :x2),[1;2]) ])
#
# plotKDE(marginal(getVertKDE(fg,:x2),[1;2]))
#
#
# @show ef2pts


























#
