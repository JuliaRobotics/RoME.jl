using IncrementalInference, TransformUtils
using Base.Test

using KernelDensityEstimate

N = 200
fg = emptyFactorGraph()


sqrtinv = [[10;0;0;0;0;0]';
[0;10;0;0;0;0]';
[0;0;10;0;0;0]';
[0;0;0;57.2958;0;0]';
[0;0;0;0;57.2958;0]';
[0;0;0;0;0;57.2958]'];

initCov = inv(sqrtinv^2)
odoCov = deepcopy(initCov)

sqrtinv = [[286.479;0]';
[0;200]']
brCov = inv(sqrtinv^2)


println("Adding PriorPose3 to graph...")
v1 = addNode!(fg,"x1",  0.1*randn(6,N),  N=N)
initPosePrior = PriorPose3(SE3([0.0;0.0;0.0], Euler(-0.00772052, 0.0, -0.0992321)), initCov)
f1  = addFactor!(fg,[v1], initPosePrior)


function addPose3Pose3(fgl::FactorGraph, n::ASCIIString, DX::SE3, cov::Array{Float64,2};
                  N::Int=100, ready::Int=1)

	prev, X, nextn = getLastPose2D(fgl)

	pts0 = projectParticles(getVal(fgl,"x1"), DX, cov)
	DXconstr = Pose3Pose3(DX, cov)
	v = addNode!(fgl, n,  pts0, N=N)
	f = addFactor!(fgl,[prev;v], DXconstr)

	return v, f
end


# Odometry pose 0 pose 1: (-0.380711, -1.02585, 1.77348; -2.19796, -0.151721, -0.0929671)
odo = SE3([-0.380711, -1.02585, 1.77348], Euler(-0.0929671, -0.151721, -2.19796) )
v, f = addPose3Pose3(fg, "x2", deepcopy(odo), odoCov, N=N)

# Odometry pose 1 pose 2: (-0.00138117, -0.0600476, 0.0204065; -0.038278, 0.0186151, 0.00606331)
odo = SE3([-0.00138117; -0.0600476; 0.0204065], Euler(0.00606331, 0.0186151, -0.038278) )
v, f = addPose3Pose3(fg, "x3", deepcopy(odo), odoCov, N=N)

# Odometry pose 2 pose 3: (-0.298704, -0.113974, 0.0244868; -0.00805163, -0.00425057, -0.0275746)
odo = SE3([-0.298704; -0.113974; 0.0244868], Euler(-0.0275746, -0.00425057, -0.00805163) )
v, f = addPose3Pose3(fg, "x4", deepcopy(odo), odoCov, N=N)

# Odometry pose 3 pose 4: (0.670492, 0.0251882, -0.0705555; 0.109113, 0.0228966, 0.039246)
odo = SE3([0.670492; 0.0251882; -0.0705555], Euler(0.039246, 0.0228966, 0.109113) )
v, f = addPose3Pose3(fg, "x5", deepcopy(odo), odoCov, N=N)

# Odometry pose 4 pose 5: (-0.0615969, 0.00849829, -0.0170747; 0.0315662, -0.0115416, -0.00605423)
odo = SE3([-0.0615969; 0.00849829; -0.0170747], Euler( -0.00605423,-0.0115416,0.0315662) )
v, f = addPose3Pose3(fg, "x6", deepcopy(odo), odoCov, N=N)


tree = prepBatchTree!(fg);
[inferOverTree!(fg, tree) for i in 1:3]; # should not be required

draw(PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
 plotKDE( getVertKDE(fg,"x6"), dimLbls=["x";"y";"z";"phi";"the";"psi"]) )


# # new install
# Pkg.add("KernelDensityEstimate")
#
# # or update
# Pkg.update()
#
# # draw this
# using KernelDensityEstimate, Gadfly
#
# p, q = kde!(randn(2,100)), kde!(randn(2,100)+2);
#
# plotKDE(p)
# plotKDE([p;q],c=["red";"blue"],levels=3)
#
# p, q = kde!(randn(4,100)), kde!(randn(4,100)+2);
# draw(PDF("test.pdf",30cm,20cm),
#         plotKDE( [p;q],c=["red";"blue"],levels=3,
# 				dimLbls=["x";"y";"z";"phi";"the";"psi"]) )









#
