# test sonar functions directly

using RoME
using Distributions

##

@testset "bearing range elevation functions" begin

##

meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

N=100
X = 0.01*randn(6,N)


# pre-emptively populate the measurements, kept separate since nlsolve calls fp(x, res) multiple times
# measurement = getSample(meas, N)

fg = initfg()
X0 = addVariable!(fg, :x0, Pose3)
initManual!(fg, :x0, X)
X1 = addVariable!(fg, :x1, Point3)
addFactor!(fg, [:x0;:x1], meas, graphinit=false)

pts = approxConv(fg, :x0x1f1, :x1)

##

@warn "still need to insert kld(..) test to ensure this is working"

p1 = kde!(pts);


println("Test back projection from ")

meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

N = 100
L = zeros(3,N);
L[1,:] .+= 3.0
L[2,:] .+= 0.65

fg = initfg()
X0 = addVariable!(fg, :x0, Pose3)
X1 = addVariable!(fg, :x1, Point3)
initManual!(fg, :x1, L)
addFactor!(fg, [:x0;:x1], meas, graphinit=false)

pts = approxConv(fg, :x0x1f1, :x0)


p2 = kde!(pts);


end

# using Gadfly
#
# axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )
# #
#
# axis = [[-2.5;2.5]';[-2.5;2.5]';[-2.5;2.5]';[-2pi;2pi]';[-2pi;2pi]';[-2pi;2pi]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p2, dimLbls=["x";"y";"z"; "roll"; "pitch"; "yaw"], axis=axis)  )


##

@testset "LinearRangeBearingElevation in solve" begin

##

# do test directly in factor graph
fg = initfg()
# @show fg.registeredModuleFunctions
N = 100
@warn "Breaks if not set to 100"

initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)


v1 = addVariable!(fg,:x1, Pose3,  N=N)
initPosePrior = PriorPose3( MvNormal(zeros(6), initCov) )
f1  = addFactor!(fg,[v1], initPosePrior)

# implemented in SensorModels.jl
meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

@time X = getVal(v1)

v2 = addVariable!(fg, :l1, Point3, N=N)
f2 = addFactor!(fg, [:x1;:l1], meas) #, threadmodel=MultiThreaded)

# L1pts = approxConv(fg, :x1l1f1, :l1)
# X1pts = approxConv(fg, :x1l1f1, :x1)

solveTree!(fg);

##

end