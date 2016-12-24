# test sonar functions directly

# addprocs(2)

using KernelDensityEstimate
using TransformUtils
using IncrementalInference
using RoME


# Functor for efficient functional programming, avoids type_inference at each call
fp! = WrapParam{reuseLBRA}(zeros(3), zeros(6), zeros(3), reuseLBRA(0))

# @everywhere begin
# (p::WrapParam{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
#     residualLRBE!(res, p.z, x, p.landmark, p.reuse)
# end

meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))


@time fp!(zeros(6), zeros(3))
fgr = FastGenericRoot{WrapParam}(6, 3, fp!)

# LinearRangeBearingElevation(Distributions.Normal{Float64}(μ=0.2, σ=0.001),Distributions.Normal{Float64}(μ=3.0, σ=0.001),Distributions.Uniform{Float64}(a=-0.25133, b=0.25133))
X, pts = 0.01*randn(6,200), zeros(3,200);
# X[1:3,:] = 0.1*randn(3,200)
@time for i in 1:200
	project!(meas, X, pts, i, fp!)
end

p1 = kde!(pts);


pts, L = 0.01*randn(6,200), zeros(3,200);
L[1,:] += 3.0
L[2,:] += 0.65


# using ProfileView
# Profile.clear()

@time for i in 1:200
  backprojectRandomized!(meas, L, pts, i, fp!)
end

# ProfileView.view()


# work on speeding up and refactoring the numericroot operation
# did not yet have the desired affect, lots of memory still being claimed
println("Compute with FastGenericRoot memory structure")
fgr = FastGenericRoot{WrapParam}(6, 3, fp!)
@time for idx in 1:200
  fp!.landmark[1:3] = L[1:3,idx]
  fp!.pose[1:6] = pts[1:6,idx]
  fp!.z[1:3] = getSample(meas)
  copy!(fgr.X, pts[1:6,idx]) #initial guess x0
  backprojectRandomized!(fgr)
	pts[1:6,idx] = fgr.Y[1:6]
end


pts, L = 0.01*randn(6,200), zeros(3,200);
L[1,:] += 3.0
L[2,:] += 0.65


# going faster
fpA! = WrapParamArray(L, pts, zeros(3), 1, reuseLBRA(0))

# @everywhere begin
# (p::WrapParamArray{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
#     residualLRBE!(res, p.z, x, p.landmark[:,p.idx], p.reuse)
# end
# Profile.clear()

println("Compute with FastGenericRoot memory structure")
fgr = FastGenericRoot{WrapParamArray}(6, 3, fpA!)

@time for idx in 1:200
  fpA!.idx = idx
  getSample!(fpA!.z,  meas)
  # fpA!.z[1:3] = getSample(meas)
  copy!(fgr.X, pts[1:6,idx]) #initial guess x0
  backprojectRandomized!( fgr )
	pts[1:6,idx] = fgr.Y[1:6]
end


# ProfileView.view()


p2 = kde!(pts);




#  visualization

#
# using Gadfly
#
# axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )
#
#
# axis = [[-2.5;2.5]';[-2.5;2.5]';[-2.5;2.5]';[-2pi;2pi]';[-2pi;2pi]';[-2pi;2pi]']
# draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
#       plotKDE( p2, dimLbls=["x";"y";"z"; "roll"; "pitch"; "yaw"], axis=axis)  )



# do test directly in factor graph
fg = initfg()
N = 200

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
for i in 1:N
	project!(meas, X, pts, i, fp!)
end
v2 = addNode!(fg,:l1,  pts, N=N)
f2 = addFactor!(fg,[v1;v2],meas)


L1pts = evalFactor2(fg, f2, fg.IDs[:l1])
X1pts = evalFactor2(fg, f2, fg.IDs[:x1])



tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree)
inferOverTree!(fg, tree)
