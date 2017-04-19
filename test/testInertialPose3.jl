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



v2 = addNode!(fg, :x2, dims=15,  N=N)

n = 1
noise = MvNormal(zeros(15),(DATA[n][3]+DATA[n][3]')*0.5 )
inerodo = InertialPose3(noise,DATA[n][1],DATA[n][2])

f2  = addFactor!(fg, [v1;v2], inerodo )


initializeNode!(fg, :x2, N=N)



using Gadfly

spy(DATA[n][3])

norm(DATA[n][3]-DATA[n][3]')



tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree)




plotKDE(fg, :x1, dims=[4;5])





writeGraphPdf(fg)
run(`evince fg.pdf`)










#
