# test packing functions


using RoME
# using IncrementalInference

using Base.Test


N = 100
fg = initfg()


initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.01])

# Some starting position
v1 = addNode!(fg, :x1, zeros(3,1), diagm([1.0;1.0;0.1]), N=N)
ipp = PriorPose2(zeros(3,1), initCov, [1.0])
f1  = addFactor!(fg,[v1], ipp)

# and a second pose
v2 = addNode!(fg, :x2, ([50.0;0.0;pi/2]')', diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;pi/2]')', odoCov, [1.0])
f2 = addFactor!(fg, [v1;v2], ppc)


println("test conversions of PriorPose2")
dd = convert(PackedPriorPose2, ipp)
upd = convert(RoME.PriorPose2, dd)

@test compare(ipp, upd)

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose2}, getData(f1))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{RoME.PriorPose2}, packeddata)

@test IncrementalInference.compare(getData(f1), unpackeddata)

packedv1data = VNDencoder(IncrementalInference.PackedVariableNodeData, getData(v1))
upv1data = VNDdecoder(IncrementalInference.VariableNodeData, packedv1data)

@test IncrementalInference.compare(getData(v1), upv1data)



println("test conversions of Pose2Pose2")

dd = convert(PackedPose2Pose2, ppc)
upd = convert(RoME.Pose2Pose2, dd)

# TODO -- fix ambiguity in compare function
@test compare(ppc, upd)

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Pose2}, getData(f2))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{RoME.Pose2Pose2}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f2), unpackeddata)





# parameters
N = 300
initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)

# start new factor graph
fg = initfg()

v1 = addNode!(fg,:x1,  0.1*randn(6,N))
ipp = PriorPose3(SE3(0), initCov)
f1  = addFactor!(fg,[v1], ipp)


println("test conversions of PriorPose3")
dd = convert(PackedPriorPose3, ipp)
upd = convert(RoME.PriorPose3, dd)

# TODO -- fix ambibuity in compare function
@test TransformUtils.compare(ipp.Zi, upd.Zi)
@test norm(ipp.Cov - upd.Cov) < 1e-8

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose3}, getData(f1))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{RoME.PriorPose3}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f1), unpackeddata)


packedv1data = VNDencoder(IncrementalInference.PackedVariableNodeData, getData(v1))
upv1data = VNDdecoder(IncrementalInference.VariableNodeData, packedv1data)

@test IncrementalInference.compare(getData(v1), upv1data)




println("test conversions of Pose3Pose3")

pp3 = Pose3Pose3( SE3([25;0;0], Quaternion(0)), odoCov)
v2, f2 = addOdoFG!(fg, pp3 )

dd = convert(PackedPose3Pose3, pp3)
upd = convert(RoME.Pose3Pose3, dd)

# TODO -- fix ambiguity in compare function
@test TransformUtils.compare(pp3.Zij, upd.Zij)
@test norm(pp3.Cov - upd.Cov) < 1e-8

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3}, getData(f2))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{RoME.Pose3Pose3}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f2), unpackeddata)




println("test conversions of Pose3Pose3NH")

odo = SE3(randn(3), convert(Quaternion, so3(0.1*randn(3)) ))
odoc = Pose3Pose3NH(odo, odoCov, [0.5;0.5])

f3 = addFactor!(fg,[:x1;:x2],odoc)


dd = convert(PackedPose3Pose3NH, odoc)
upd = convert(RoME.Pose3Pose3NH, dd)

@test TransformUtils.compare(odoc.Zij, upd.Zij)
@test norm(odoc.Cov - upd.Cov) < 1e-8
@test norm(odoc.ValidHypot.p - upd.ValidHypot.p) < 1e-8


packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3NH}, getData(f3))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{RoME.Pose3Pose3NH}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f3), unpackeddata)





















#
