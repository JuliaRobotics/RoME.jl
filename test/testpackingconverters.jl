# test packing functions


using RoME
using Distributions
using TransformUtils
using KernelDensityEstimate

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

@test RoME.compare(ipp, upd)

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose2}, getData(f1))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.PriorPose2}}, packeddata)

@test IncrementalInference.compare(getData(f1), unpackeddata)

packedv1data = VNDencoder(IncrementalInference.PackedVariableNodeData, getData(v1))
upv1data = VNDdecoder(IncrementalInference.VariableNodeData, packedv1data)

@test IncrementalInference.compare(getData(v1), upv1data)



println("test conversions of Pose2Pose2")

dd = convert(PackedPose2Pose2, ppc)
upd = convert(RoME.Pose2Pose2, dd)

# TODO -- fix ambiguity in compare function
@test RoME.compare(ppc, upd)

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Pose2}, getData(f2))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose2Pose2}}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f2), unpackeddata)


println("test conversions of Pose2DPoint2DBearingRange")


# and a second pose
v3 = addNode!(fg, :l1, ([50.0;50.0]')', diagm([1.0;1.0]), N=N)
# ppc = Pose2DPoint2DBearingRange([50.0;0.0;pi/2], 0.01*eye(2), [1.0])
ppbr = Pose2DPoint2DBearingRange{Normal, Normal}(
              Normal(0.0, 0.005 ),
              Normal(50, 0.5) )
f3 = addFactor!(fg, [:x2;:l1], ppbr)

dd = convert(PackedPose2DPoint2DBearingRange, ppbr)
upd = convert(RoME.Pose2DPoint2DBearingRange{Normal,Normal}, dd)


packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2DPoint2DBearingRange}, getData(f3))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose2DPoint2DBearingRange{Normal,Normal}}}, packeddata)

@test ppbr.bearing.μ == unpackeddata.fnc.usrfnc!.bearing.μ
@test ppbr.bearing.σ == unpackeddata.fnc.usrfnc!.bearing.σ

@test ppbr.range.μ == unpackeddata.fnc.usrfnc!.range.μ
@test ppbr.range.σ == unpackeddata.fnc.usrfnc!.range.σ


println("test conversions of PriorPose3")

# parameters
N = 300
initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)

# start new factor graph
fg = initfg()

v1 = addNode!(fg,:x1,  0.1*randn(6,N))
ipp = PriorPose3(MvNormal(zeros(6), initCov) )
f1  = addFactor!(fg,[v1], ipp)

dd = convert(PackedPriorPose3, ipp)
upd = convert(RoME.PriorPose3, dd)

# @test TransformUtils.compare(ipp.Zi, upd.Zi)
@test norm(ipp.Zi.μ - upd.Zi.μ) < 1e-10
@test norm(ipp.Zi.Σ.mat - upd.Zi.Σ.mat) < 1e-8

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose3}, getData(f1))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.PriorPose3}}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f1), unpackeddata)

packedv1data = VNDencoder(IncrementalInference.PackedVariableNodeData, getData(v1))
upv1data = VNDdecoder(IncrementalInference.VariableNodeData, packedv1data)

@test IncrementalInference.compare(getData(v1), upv1data)




println("test conversions of Pose3Pose3")

pp3 = Pose3Pose3( MvNormal([25.0;0;0;0;0;0], odoCov) )
v2, f2 = addOdoFG!(fg, pp3 )

dd = convert(PackedPose3Pose3, pp3)
upd = convert(RoME.Pose3Pose3, dd)


@test norm(pp3.Zij.μ - upd.Zij.μ) < 1e-10
@test norm(pp3.Zij.Σ.mat - upd.Zij.Σ.mat) < 1e-8

packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3}, getData(f2))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose3Pose3}}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f2), unpackeddata)




println("test conversions of Pose3Pose3NH")

odo = SE3(randn(3), convert(Quaternion, so3(0.1*randn(3)) ))
odoc = Pose3Pose3NH( MvNormal(veeEuler(odo),odoCov), [0.5;0.5])

f3 = addFactor!(fg,[:x1;:x2],odoc)


dd = convert(PackedPose3Pose3NH, odoc)
upd = convert(RoME.Pose3Pose3NH, dd)

@test norm(odoc.Zij.μ - upd.Zij.μ) < 1e-8
@test norm(odoc.Zij.Σ.mat - upd.Zij.Σ.mat) < 1e-8
@test norm(odoc.nullhypothesis.p - upd.nullhypothesis.p) < 1e-8


packeddata = FNDencode(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3NH}, getData(f3))
unpackeddata = FNDdecode(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose3Pose3NH}}, packeddata)

# TODO -- fix ambibuity in compare function
@test IncrementalInference.compare(getData(f3), unpackeddata)




println("test conversions of PartialPriorRollPitchZ")

prpz = PartialPriorRollPitchZ(MvNormal([0.0;0.5],0.1*eye(2)),Normal(3.0,0.5))

pprpz = convert(PackedPartialPriorRollPitchZ, prpz)
unp = convert(PartialPriorRollPitchZ, pprpz)

@test RoME.compare(prpz, unp)





println("test conversions of PartialPose3XYYaw")


xyy = PartialPose3XYYaw(MvNormal([1.0;2.0;0.5],0.1*eye(3)))

pxyy = convert(PackedPartialPose3XYYaw, xyy)
unp = convert(PartialPose3XYYaw, pxyy)

@test RoME.compare(xyy, unp)




println("test conversions of PartialPose3XYYawNH")


xyy = PartialPose3XYYawNH(MvNormal([1.0;2.0;0.5],0.1*eye(3)), [0.6;0.4])

pxyy = convert(PackedPartialPose3XYYawNH, xyy)
unp = convert(PartialPose3XYYawNH, pxyy)

@test RoME.compare(xyy, unp)



println("test PriorPoint2DensityNH")

prpt2 = PriorPoint2DensityNH(kde!(randn(2,100)),[0.25;0.75]  )

pprpt2 = convert(PackedPriorPoint2DensityNH, prpt2)
uprpt2 = convert(PriorPoint2DensityNH, pprpt2)

@test norm(getPoints(prpt2.belief)-getPoints(uprpt2.belief)) < 1e-8

@test norm(prpt2.nullhypothesis.p-uprpt2.nullhypothesis.p) < 1e-8







#
