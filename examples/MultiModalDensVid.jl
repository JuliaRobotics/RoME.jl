using RoME, TransformUtils, KernelDensityEstimate
using IncrementalInference, CloudGraphs, Neo4j


include("BlandAuthDB.jl")

configuration = CloudGraphs.CloudGraphConfiguration(dbaddress, 7474, dbusr, dbpwd, mongoaddress, 27017, false, "", "");
cloudGraph = connect(configuration, IncrementalInference.getfnctype);
# register types of interest in CloudGraphs
registerGeneralVariableTypes!(cloudGraph)
IncrementalInference.setCloudDataLayerAPI!()


N = 200
fg = initfg()
fg.sessionname="SESSANIM"
fg.cg = cloudGraph


initCov = Matrix(Diagonal([0.03;0.03;0.001].^2))
odoCov = Matrix(Diagonal([3.0;3.0;0.4].^2))

x1 = initFactorGraph!(fg, ready=0, N=N)
println("Starting $(fg.sessionname) at $(x1), at $(getKDEMax(getVertKDE(fg,x1)))")

r = 3.0/12.0*pi
cDX = [10.0;0;r]
for i in 2:13
	addOdoFG!(fg, "x$(i)", cDX, odoCov, ready=0, N=N)
end


# insert loop closure
v1a = getVert(fg, "x1")
v1b = getVert(fg, "x13")

DX = zeros(3,1)
pp = Pose2Pose2(DX, 0.1*Matrix{Float64}(LinearAlgebra.I, 3,3), [1.0]) #[prev;v],
f = addFactor!(fg, [v1a;v1b], pp, ready=0)

v1c = getVert(fg, "x7")

DX = zeros(3,1)
DX[2] = 38.63
DX[3] = pi
pp = Pose2Pose2(DX, 0.1*Matrix{Float64}(LinearAlgebra.I, 3,3), [1.0]) #[prev;v],
f = addFactor!(fg, [v1a;v1c], pp, ready=0)


v2a = getVert(fg, "x4")
v2b = getVert(fg, "x10")
pp = Pose2Pose2(DX, 0.1*Matrix{Float64}(LinearAlgebra.I, 3,3), [0.1;0.9]) #[prev;v],
f = addFactor!(fg, [v2a;v2b], pp, ready=0)


setDBAllReady!(fg)




#
