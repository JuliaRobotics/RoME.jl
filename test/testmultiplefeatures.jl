#  test MultipleFeatures2D
using RoME
# , Distributions
# using IncrementalInference   # for evalFactor2
using Test

using Optim

# using Gadfly

global x1 = [0.0;0;0]
global x2 = [1.0;-1.0;pi/2]

global l1 = [0.5;0.5]
global l2 = [1.0;0.0]
global l3 = [1.5;-0.5]
global l3b = [1.5;0.0]
global L = [l1';l2';l3']

# plot(x=L[:,1],y=L[:,2],Geom.point)

global the11 = atan(l1[2],l1[1])
global the12 = atan(l2[2],l2[1])
global the13 = atan(l3[2],l3[1])

global Xi = SE2(x1)
global Xj = SE2(x2)
global Δxij = se2vee(Xi\Xj)

global tl1 = se2vee(Xj \ SE2([l1;0.0]))[1:2]
global the21 = atan(tl1[2],tl1[1])

global tl2 = se2vee(Xj \ SE2([l2;0.0]))[1:2]
global the22 = atan(tl2[2],tl2[1])

global tl3 = se2vee(Xj \ SE2([l3;0.0]))[1:2]
global the23 = atan(tl3[2],tl3[1])

global tl3b = se2vee(Xj \ SE2([l3b;0.0]))[1:2]
global the23b = atan(tl3b[2],tl3b[1])

global xir1 = Normal(the11,1e-4)
global xir2 = Normal(the12,1e-4)
global xir3 = Normal(the13,1e-4)
global xjr1 = Normal(the21,1e-4)
global xjr2 = Normal(the22,1e-4)
global xjr3 = Normal(the23,1e-4)
global xjr3b = Normal(the23b,1e-4)

# Not using second feature
global mm2 = MultipleFeatures2D(
  xir1,
  xir2,
  xir3,
  xjr1,
  xjr2,
  xjr3,
  xjr3b,
  Categorical([1.0;0.0]),
  SE2(zeros(3))  )

global res = randn(3)
global idx = 1
global meas = getSample(mm2)


global wAbi = reshape(x1,3,1)
global wAbj = reshape(x2,3,1)
global wAo1 = reshape(l1,2,1)
global wAo2 = reshape(l2,2,1)
global wAo3 = reshape(l3,2,1)


@testset "test MultipleFeatures unit vector functions..." begin

global res, rhat, resid = zeros(3),zeros(2), zeros(2)

for i in 1:3
  global r1, rhathat1, α1 = getUvecScaleFeature2D(Xi, SE2(zeros(3)), SE2([L[i,:];0.0]))
  rhat[1] = cos(meas[1][i,idx,1])
  rhat[2] = sin(meas[1][i,idx,1])
  # @show rhat, rhathat1
  resid[1:2] = rhat-rhathat1
  res[1] += dot(resid, resid)
end

@test 0.0 <= res[1] < 1e-5

end


@testset "test residual function..." begin

global res = zeros(3)

wAbi[:] .= 0.0
# @show wAbi

mm2(res,
  nothing,
  idx,
  meas,
  wAbi,
  wAbj,
  wAo1,
  wAo2,
  wAo3  )

# @show res

function gg2(x::Float64,y::Float64)
  wAbi[2] = x
  wAbi[3] = y

  global res = zeros(3)

  mm2(res, nothing, idx, meas,
    wAbi,
    wAbj,
    wAo1,
    wAo2,
    wAo3  )
  return sum(res)[1]
end
global XX = range(-0.5,stop=0.5,length=200)
global YY = range(-pi/4,stop=pi/4,length=100)
gg2(0.0,0.0)
# Gadfly.plot(z=gg2,x=XX,y=YY, Geom.contour)

end


function minmickey(x::Array{Float64,1})
  wAbi[1:3,1] = x[1:3]
  res = zeros(3)
  mm2(res, nothing, idx, meas,
    wAbi,
    wAbj,
    wAo1,
    wAo2,
    wAo3  )

  return sum(res[1:2])
end


# ar99
@testset "test as minimization problem to pose..." begin


gg = (x,y) -> minmickey([x;y;0.0])
# gg = (x,y) -> x*y
# plot(z=gg, x=range(-2,stop=2,length=100),y=range(-2,stop=2,length=100), Geom.contour)


gg(0.0,0.1)

@time r = optimize(minmickey, randn(3))
@show r.minimizer

@test norm(res) < 1e-3

end


function minmickey(x::Array{Float64,1})
  wAo1[1:2,1] = x[1:2]
  res = zeros(3)
  mm2(res, nothing, idx, meas,
    wAbi,
    wAbj,
    wAo1,
    wAo2,
    wAo3  )

  return sum(res[1:2])
end


@testset "test as minimization problem to landmark..." begin

wAbi[:] .= 0.0

gg = (x,y) -> minmickey([x;y])
#plot(z=gg, x=range(1,stop=2,length=100),y=range(-1,stop=1,length=100), Geom.contour(levels=50))

# now build in factor graph form for further testing

global N=50
global fg = initfg()

global initCov = 0.05*Matrix{Float64}(LinearAlgebra.I, 3,3)
initCov[3,3] = 0.001
global odoCov = deepcopy(initCov)

global ipp = PriorPose2(MvNormal(zeros(3), initCov^2))

global v1 = addNode!(fg, :x1, Pose2, N=N)

global f1  = addFactor!(fg,[:x1], ipp)


global v2 = addNode!(fg, :x2, Pose2, N=N)

global px2 = zeros(3)
px2[1:3] = [1.1;-1.0;pi/2]
global ipp2 = PriorPose2(MvNormal(px2, initCov^2))
global f1  = addFactor!(fg,[:x2], ipp2)

global vl1 = addNode!(fg, :l1, Point2, N=N)
global vl2 = addNode!(fg, :l2, Point2, N=N)
global vl3 = addNode!(fg, :l3, Point2, N=N)

# # why is explicit call to autoinit required here? should not be necessary
# doautoinit!(fg, :x1)
# doautoinit!(fg, :x2)

global f2 = addFactor!(fg, [v1;v2;vl1;vl2;vl3], mm2, threadmodel=SingleThreaded )

# getVal(fg, :l3)
#
# writeGraphPdf(fg)

global data = getData(f2)
# fieldnames(data.fnc)
@test data.fnc.threadmodel == SingleThreaded

global ef2pts = approxConv(fg, :x1x2l1l2l3f1, :l2)
# evalFactor2(fg, f2, fg.IDs[:l2])


ensureAllInitialized!(fg)
global tree = wipeBuildNewTree!(fg)
# spyCliqMat(tree.cliques[1])

[inferOverTreeR!(fg,tree, N=N, dbg=true) for i in 1:1]



end

# lets see what is happening during MCMC runs
# plotMCMC(tree, :l2, show=true, levels=3) # make true to see pictures, false for testing

# spyCliqMat(tree, :x2)

#
# pts, = getSample(ipp2, 100);
#
# plotKDE( kde!(pts[1:2,:]) )
#
# @show ipp2.Cov
# pts = rand(MvNormal(ipp2.Zi[:,1],ipp2.Cov),N);

# using KernelDensityEstimate
# plotKDE(kde!(ef2pts[1,:]))

# pl = plotKDE([kde!(randn(2,100));kde!(2+randn(2,100))],levels=1,c=["red";"green"],legend=["1";"2"])

# arrp = [getVertKDE(fg, :l1);getVertKDE(fg, :l2);getVertKDE(fg, :l3)]
# arrp1 = deepcopy(arrp)
# plotKDE([arrp;arrp1],levels=1,c=["red";"red";"red";"blue";"blue";"blue"])

# arrp=[marginal(getVertKDE(fg, :x1),[1;2]);marginal(getVertKDE(fg, :x2),[1;2]) ]
# plotKDE(arrp,levels=1)

# plotKDE(marginal(getVertKDE(fg,:x1),[3]))

# pp,arr,parts = localProduct(fg,:x2)
# cc = plotKDE([pp;arr],c=["red";"black";"blue"],levels=1);

# cc = drawLbl(fg,:l1);
# Gadfly.draw(PDF("/home/dehann/Desktop/test.pdf",20cm,15cm),cc)





@testset "test ambiguous bi-modal multifeature constraint operation" begin

global mm2 = MultipleFeatures2D(
  xir1,
  xir2,
  xir3,
  xjr1,
  xjr2,
  xjr3,
  xjr3b,
  Categorical([0.5;0.5]),
  SE2(zeros(3))  )



end












#
