# test sonar functions directly

# addprocs(2)

using KernelDensityEstimate
using TransformUtils
using IncrementalInference
using RoME


# @everywhere begin
# function (p::LinearRangeBearingElevation)(
#       res::Vector{Float64},
#       idx::Int,
#       meas::Array{Float64,2},
#       pose::Array{Float64,2},
#       landm::Array{Float64,2}  )
#   #
#   println("LinearRangeBearingElevation is happening")
#   residualLRBE!(res, p.meas[:,idx], pose[:,idx], landm[:,idx], p.reuse)
#   nothing
# end

meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

# Functor for efficient functional programming, avoids type_inference at each call
# fp! = WrapParam{reuseLBRA}(zeros(3), zeros(6), zeros(3), reuseLBRA(0))

N = 200
X, pts = 0.01*randn(6,N), zeros(3,N);
t = Array{Array{Float64,2},1}()
push!(t,X)
push!(t,pts)

fp! = GenericWrapParam{LinearRangeBearingElevation}(meas, t, 2, 1, zeros(0,1), RoME.getSample)

# pre-emptively populate the measurements, kept separate since nlsolve calls fp(x, res) multiple times
fp!.measurement = fp!.samplerfnc(fp!.usrfnc!, N)
# fp!(x, res)
@time fp!(zeros(3), zeros(3))

@show zDim = size(fp!.measurement,1)
fr = FastRootGenericWrapParam{LinearRangeBearingElevation}(fp!.params[fp!.varidx], zDim, fp!)

@time for fp!.particleidx in 1:N
  numericRootGenericRandomizedFnc!( fr )
end
# LinearRangeBearingElevation(Distributions.Normal{Float64}(μ=0.2, σ=0.001),Distributions.Normal{Float64}(μ=3.0, σ=0.001),Distributions.Uniform{Float64}(a=-0.25133, b=0.25133))
# X[1:3,:] = 0.1*randn(3,200)
# @time for i in 1:200
# 	project!(meas, X, pts, i, fp!)
# end

warn("still need to insert kld(..) test to ensure this is working")

p1 = kde!(pts);


println("Test back projection from ")

meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))

N = 200
pts, L = 0.01*randn(6,N), zeros(3,N);
L[1,:] += 3.0
L[2,:] += 0.65
t = Array{Array{Float64,2},1}()
push!(t,pts)
push!(t,L)

fp! = GenericWrapParam{LinearRangeBearingElevation}(meas, t, 1, 1, zeros(0,1), RoME.getSample)

# pre-emptively populate the measurements, kept separate since nlsolve calls fp(x, res) multiple times
fp!.measurement = fp!.samplerfnc(fp!.usrfnc!, N)
# fp!(x, res)
@time fp!(zeros(6), zeros(3))

@show zDim = size(fp!.measurement,1)
fr = FastRootGenericWrapParam{LinearRangeBearingElevation}(fp!.params[fp!.varidx], zDim, fp!)

@time for fp!.particleidx in 1:N
  numericRootGenericRandomizedFnc!( fr )
end


p2 = kde!(pts);



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
fg = initfg()
# @show fg.registeredModuleFunctions
N = 100
warn("Breaks if not set to 100")

initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)

println("Adding PriorPose3 to graph...")
X, pts = 0.01*randn(6,N), zeros(3,N);

v1 = addNode!(fg,:x1,  X,  N=N)
initPosePrior = PriorPose3(SE3(0), initCov)
f1  = addFactor!(fg,[v1], initPosePrior)


println("Adding LinearRangeBearingElevation to graph...")
meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))
# for i in 1:N
# 	project!(meas, X, pts, i, fp!)
# end
# implemented in SensorModels
@time X = getVal(v1)
@time pts = X + meas
p1 = kde!(pts); # visual checking

v2 = addNode!(fg, :l1, pts, N=N)
f2 = addFactor!(fg,[v1;v2],meas)


L1pts = evalFactor2(fg, f2, fg.IDs[:l1])
X1pts = evalFactor2(fg, f2, fg.IDs[:x1])


tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree)
inferOverTree!(fg, tree)




# using ProfileView
# Profile.clear()

# @time for i in 1:200
#   backprojectRandomized!(meas, L, pts, i, fp!)
# end
#
# # ProfileView.view()
#
#
# # work on speeding up and refactoring the numericroot operation
# # did not yet have the desired affect, lots of memory still being claimed
# println("Compute with FastGenericRoot memory structure")
# fgr = FastGenericRoot{WrapParam}(6, 3, fp!)
# @time for idx in 1:200
#   fp!.landmark[1:3] = L[1:3,idx]
#   fp!.pose[1:6] = pts[1:6,idx]
#   fp!.z[1:3] = getSample(meas)
#   copy!(fgr.X, pts[1:6,idx]) #initial guess x0
#   numericRootGenericRandomizedFnc!( fgr )
#   # backprojectRandomized!(fgr)
# 	pts[1:6,idx] = fgr.Y[1:6]
# end
#
#
# pts, L = 0.01*randn(6,200), zeros(3,200);
# L[1,:] += 3.0
# L[2,:] += 0.65
#
#
# # going faster
# fpA! = WrapParamArray(L, pts, zeros(3), 1, reuseLBRA(0))
#
# # @everywhere begin
# # (p::WrapParamArray{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
# #     residualLRBE!(res, p.z, x, p.landmark[:,p.idx], p.reuse)
# # end
# # Profile.clear()
#
# println("Compute with FastGenericRoot memory structure")
# fgr = FastGenericRoot{WrapParamArray}(6, 3, fpA!)
#
# @time for idx in 1:200
#   fpA!.idx = idx
#   getSample!(fpA!.z,  meas)
#   # fpA!.z[1:3] = getSample(meas)
#   copy!(fgr.X, pts[1:6,idx]) #initial guess x0
#   numericRootGenericRandomizedFnc!( fgr )
#   # backprojectRandomized!( fgr )
# 	pts[1:6,idx] = fgr.Y[1:6]
# end
#

# ProfileView.view()



#  visualization
