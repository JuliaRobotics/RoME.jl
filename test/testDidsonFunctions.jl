# test sonar functions directly

# addprocs(2)

# using KernelDensityEstimate
# using TransformUtils
# using IncrementalInference
using RoME
using Distributions



global meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))


global N = 100
global X, pts = 0.01*randn(6,N), zeros(3,N);
global t = Array{Array{Float64,2},1}()
push!(t,X)
push!(t,pts)


# pre-emptively populate the measurements, kept separate since nlsolve calls fp(x, res) multiple times
global measurement = getSample(meas, N)

@show zDim = size(measurement[1],1)

global ccw = CommonConvWrapper(meas, t[2], zDim, t, measurement=measurement, varidx=2)

# TODO remove
ccw.measurement = measurement
ccw.cpt[1].res = zeros(1)

@time ccw(zeros(3), zeros(3))
@time ccw(zeros(3))



@time for n in 1:N
  ccw.cpt[Threads.threadid()].particleidx = n
  numericRootGenericRandomizedFnc!( ccw ) #fr
end



@warn "still need to insert kld(..) test to ensure this is working"

global p1 = kde!(pts);


println("Test back projection from ")

global meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

global N = 200
global pts, L = 0.01*randn(6,N), zeros(3,N);
L[1,:] .+= 3.0
L[2,:] .+= 0.65
global t = Array{Array{Float64,2},1}()
push!(t,pts)
push!(t,L)

global measurement = getSample(meas, N)
global zDim = size(measurement,1)
global ccw = CommonConvWrapper(meas, t[1], zDim, t, varidx=1, measurement=measurement)

# TODO remove
ccw.measurement = measurement
ccw.cpt[1].res = zeros(1)

# pre-emptively populate the measurements, kept separate since nlsolve calls fp(x, res) multiple times
@time ccw(zeros(3), zeros(6))
@time ccw(zeros(6))


@time for n in 1:N
  ccw.cpt[Threads.nthreads()].particleidx = n
  numericRootGenericRandomizedFnc!( ccw )
end


global p2 = kde!(pts);



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



# do test directly in factor graph
global fg = initfg()
# @show fg.registeredModuleFunctions
global N = 100
@warn "Breaks if not set to 100"

global initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.01 for i in 4:6];
global odoCov = deepcopy(initCov)

println("Adding PriorPose3 to graph...")
# X, pts = 0.01*randn(6,N), zeros(3,N);

global v1 = addNode!(fg,:x1, Pose3,  N=N)
global initPosePrior = PriorPose3( MvNormal(zeros(6), initCov) )
global f1  = addFactor!(fg,[v1], initPosePrior)


println("Adding LinearRangeBearingElevation to graph...")
# implemented in SensorModels
global meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

@time X = getVal(v1)
# @time pts = approxConvBinary(X, meas, 3)  # TODO add back when IIF v0.3.8+ is available
# p1 = kde!(pts); # visual checking

global v2 = addNode!(fg, :l1, Point3, N=N)
global f2 = addFactor!(fg, [:x1;:l1], meas) #, threadmodel=MultiThreaded)

# ensureAllInitialized!(fg)
# getVal(fg, :x1)


global L1pts = approxConv(fg, :x1l1f1, :l1)


global data = getData(f2)
# data


global X1pts = approxConv(fg, :x1l1f1, :x1)

# isInitialized(fg, :l1)
# getVal(fg, :l1)

ensureAllInitialized!(fg)
global tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree, N=N)
# inferOverTree!(fg, tree, N=N)
