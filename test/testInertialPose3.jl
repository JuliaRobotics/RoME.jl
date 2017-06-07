# test InertialPose3

using RoME, Distributions
using JLD


@load joinpath(dirname(@__FILE__),"preintstationarydata.jld")
# DATA[1]

N = 100
fg = initfg()


initCov = 0.001*eye(15)
[initCov[i,i] = 0.0001^2 for i in 4:6];
[initCov[i,i] = 0.0002^2 for i in 10:15];
odoCov = deepcopy(initCov)



println("Adding PriorInertialPose3 to graph...")
v1 = addNode!(fg, :x1,  0.1*randn(15,N),  N=N)
initPosePrior = PriorInertialPose3( MvNormal( zeros(15), initCov) )
f1  = addFactor!(fg, [v1], initPosePrior)


n = 1
v2 = addNode!(fg, :x2, dims=15,  N=N)
noise = MvNormal(zeros(15),(DATA[n][3]+DATA[n][3]')*0.5 )
inerodo = InertialPose3(noise,DATA[n][1],DATA[n][2])
f2  = addFactor!(fg, [v1;v2], inerodo )


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


res = zeros(15)
idx = 1
# pos, so3, vel, bw, ba
meas = ([0;0;9.81/2; 0;0;0; 0;0;9.81; 0;0;0; 0;0;0]')'
wIPi = zeros(15,1)
wIPj = zeros(15,1)

inerodo(res, idx, (meas,), wIPi, wIPj)
# function (ip3::InertialPose3)(
#         res::Vector{Float64},
#         idx::Int,
#         meas::Tuple,
#         wIPi::Array{Float64,2},
#         wIPj::Array{Float64,2}  )

using Optim

ggo = (x) -> inerodo(res, idx, (meas,), wIPi, (x')')

ggos = (x) -> ggo([x[1:5]...,0.0,x[6:end]...])

ret = optimize(ggo, zeros(15))



using Gadfly

obj = (x) -> ggo([0,0.0,0, 0,0,x, 0,0,0, 0,0,0, 0,0,0])
plot(obj, -2, 2)

ran = linspace(-2,2,100)
ggoxy = (x,y) -> ggo([0,0,0.0, 0,0,x, y,0,0, 0,0,0, 0,0,0])
plot(z=ggoxy, x=ran, y=ran, Geom.contour)


@show res
using NLsolve

gg = (x, res) -> inerodo(res, idx, (meas,), wIPi, (x')')

ret = nlsolve(gg, wIPj[:])

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



tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree)




plotKDE(fg, :x1, dims=[4;5])





writeGraphPdf(fg)
run(`evince fg.pdf`)










#
