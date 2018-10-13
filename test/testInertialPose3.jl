# test InertialPose3

using RoME
# , Distributions
using JLD


@load joinpath(dirname(@__FILE__),"preintstationarydata.jld")
# DATA[1]

global N = 100
global fg = initfg()


global initCov = 0.001*Matrix{Float64}(LinearAlgebra.I, 15,15)
[initCov[i,i] = 0.0001^2 for i in 4:6];
[initCov[i,i] = 0.0002^2 for i in 10:15];
global odoCov = deepcopy(initCov)



println("Adding PriorInertialPose3 to graph...")
global v1 = addNode!(fg, :x1, InertialPose3, N=N) #0.1*randn(15,N)
global initPosePrior = PriorInertialPose3( MvNormal( zeros(15), initCov) )
global f1  = addFactor!(fg, [v1], initPosePrior)


global n = 1
global v2 = addNode!(fg, :x2, InertialPose3, dims=15,  N=N)
global noise = MvNormal(zeros(15),(DATA[n][3]+DATA[n][3]')*0.5 )
global inerodo = InertialPose3(noise,DATA[n][1],DATA[n][2])
global f2  = addFactor!(fg, [v1;v2], inerodo )


initializeNode!(fg, :x2, N=N)

# drawDensityMesh(fg, :x2)

plotKDE(fg, :x2, dims=[1;2]  , title="x,y")
plotKDE(fg, :x2, dims=[2;3]  , title="y,z")
plotKDE(fg, :x2, dims=[6;3]  , title="ψ,z")
plotKDE(fg, :x2, dims=[4;5]  , title=",")
plotKDE(fg, :x2, dims=[7;8]  , title=",")
plotKDE(fg, :x2, dims=[9;10] , title=",")
plotKDE(fg, :x2, dims=[11;12], title=",")
plotKDE(fg, :x2, dims=[13;14], title=",")
plotKDE(fg, :x2, dims=[15]   , title=",")


DATA[n][3][1,1]


global res = zeros(15)
global idx = 1
# pos, so3, vel, bw, ba
global meas = vectoarr2([0;0;9.81/2; 0;0;0; 0;0;9.81; 0;0;0; 0;0;0])
global wIPi = zeros(15,1)
global wIPj = zeros(15,1)

inerodo(res, nothing, idx, (meas,), wIPi, wIPj)

using Optim

ggo = (x) -> inerodo(res, nothing, idx, (meas,), wIPi, vectoarr2(x))

ggos = (x) -> ggo([x[1:5]...,0.0,x[6:end]...])

global ret = optimize(ggo, zeros(15))



using Gadfly

global obj = (x) -> ggo([0,0.0,0, 0,0,x, 0,0,0, 0,0,0, 0,0,0])
plot(obj, -2, 2)

global ran = range(-2,stop=2,length=100)
ggoxy = (x,y) -> ggo([0,0,0.0, 0,0,x, y,0,0, 0,0,0, 0,0,0])
plot(z=ggoxy, x=ran, y=ran, Geom.contour)


@show res
using NLsolve

gg = (res, x) -> inerodo(res, nothing, idx, (meas,), wIPi, vectoarr2(x))

global ret = nlsolve(gg, wIPj[:])

@show ret.zero




using KernelDensityEstimate
p = kde!( getSample(inerodo, 100)[1] )
plotKDE(p, dims=[1;2]  )
plotKDE(p, dims=[6;3]  )
plotKDE(p, dims=[4;5]  )
plotKDE(p, dims=[7;8]  )
plotKDE(p, dims=[9;10] )
plotKDE(p, dims=[11;12])
plotKDE(p, dims=[13;14])
plotKDE(p, dims=[15]   )




# Juno.breakpoint("/home/dehann/.julia/v0.5/IncrementalInference/src/ApproxConv.jl", 22)

using Gadfly

spy(DATA[n][3])

norm(DATA[n][3]-DATA[n][3]')


ensureAllInitialized!(fg)
global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)




plotKDE(fg, :x1, dims=[4;5])





writeGraphPdf(fg)
run(`evince fg.pdf`)
Base.rm("fg.pdf")









#
