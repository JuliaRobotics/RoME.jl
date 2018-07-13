# test packing functions


using Distributions
using KernelDensityEstimate
using TransformUtils
using RoME

using Base.Test




@testset "test PriorPoint2" begin

  prpt2 = PriorPoint2( MvNormal([0.25;0.75], diagm([1.0;2.0]))  )

  pprpt2 = convert(PackedPriorPoint2, prpt2)
  uprpt2 = convert(PriorPoint2, pprpt2)

  @test norm(prpt2.Z.μ - uprpt2.Z.μ) < 1e-8

  @test norm(prpt2.Z.Σ.mat - uprpt2.Z.Σ.mat) < 1e-8

  # test backwards compatibility, TODO remove
  prpt2 = PriorPoint2D([0.25;0.75], diagm([1.0;2.0]), [1.0;]  )

end




N = 100
fg = initfg()


initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.01])

# Some starting position
v1 = addNode!(fg, :x1, Pose2, N=N) # zeros(3,1), diagm([1.0;1.0;0.1])
ipp = PriorPose2(MvNormal(zeros(3), initCov))
f1  = addFactor!(fg,[v1], ipp)

# and a second pose
v2 = addNode!(fg, :x2, Pose2, N=N) # vectoarr2([50.0;0.0;pi/2]), diagm([1.0;1.0;0.05])
ppc = Pose2Pose2([50.0;0.0;pi/2], odoCov)
f2 = addFactor!(fg, [:x1;:x2], ppc)



@testset "test conversions of PriorPose2" begin
    dd = convert(PackedPriorPose2, ipp)
    upd = convert(RoME.PriorPose2, dd)

    @test RoME.compare(ipp, upd) # temp use of RoME.compare

    packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose2}, getData(f1))
    unpackeddata = convert(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.PriorPose2}}, packeddata)

    # getData(f1)
    # unpackeddata

    @test RoME.compare(getData(f1), unpackeddata) # temp use of RoME.compare

    packedv1data = convert(IncrementalInference.PackedVariableNodeData, getData(v1))
    upv1data = convert(IncrementalInference.VariableNodeData, packedv1data)

    @test IncrementalInference.compare(getData(v1), upv1data)
end


@testset "test conversions of Pose2Pose2" begin

    dd = convert(PackedPose2Pose2, ppc)
    upd = convert(RoME.Pose2Pose2, dd)

    # TODO -- fix ambiguity in compare function
    @test RoME.compare(ppc, upd)

    packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Pose2}, getData(f2))
    unpackeddata = convert(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose2Pose2}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f2), unpackeddata)
end


@testset "test conversions of Pose2Point2BearingRange" begin
    # and a second pose
    v3 = addNode!(fg, :l1, Point2, N=N) # vectoarr2([50.0,50.0]), diagm([1.0;1.0])
    # ppc = Pose2Point2BearingRange([50.0;0.0;pi/2], 0.01*eye(2), [1.0])
    ppbr = Pose2Point2BearingRange(
                  Normal(0.0, 0.005 ),
                  Normal(50, 0.5) )
    f3 = addFactor!(fg, [:x2;:l1], ppbr)

    dd = convert(PackedPose2Point2BearingRange, ppbr)
    upd = convert(
            RoME.Pose2Point2BearingRange,
            dd
          )


    packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Point2BearingRange}, getData(f3))
    unpackeddata = convert(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose2Point2BearingRange}}, packeddata) # IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose2Point2BearingRange{Normal{Float64},Normal{Float64}}}}

    @test ppbr.bearing.μ == unpackeddata.fnc.usrfnc!.bearing.μ
    @test ppbr.bearing.σ == unpackeddata.fnc.usrfnc!.bearing.σ

    @test ppbr.range.μ == unpackeddata.fnc.usrfnc!.range.μ
    @test ppbr.range.σ == unpackeddata.fnc.usrfnc!.range.σ
end

# parameters
N = 300
initCov = eye(6)
[initCov[i,i] = 0.01 for i in 4:6];
odoCov = deepcopy(initCov)

# start new factor graph
fg = initfg()

v1 = addNode!(fg, :x1, Pose3, N=N) #  0.1*randn(6,N)
ipp = PriorPose3(MvNormal(zeros(6), initCov) )
f1  = addFactor!(fg,[v1], ipp)


@testset "test conversions of PriorPose3" begin

    dd = convert(PackedPriorPose3, ipp)
    upd = convert(RoME.PriorPose3, dd)

    # @test TransformUtils.compare(ipp.Zi, upd.Zi)
    @test norm(ipp.Zi.μ - upd.Zi.μ) < 1e-10
    @test norm(ipp.Zi.Σ.mat - upd.Zi.Σ.mat) < 1e-8

    packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose3}, getData(f1))
    unpackeddata = convert(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.PriorPose3}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f1), unpackeddata)

    packedv1data = convert(IncrementalInference.PackedVariableNodeData, getData(v1))
    upv1data = convert(IncrementalInference.VariableNodeData, packedv1data)

    @test IncrementalInference.compare(getData(v1), upv1data)
end



pp3 = Pose3Pose3( MvNormal([25.0;0;0;0;0;0], odoCov) )
v2, f2 = addOdoFG!(fg, pp3 )

@testset "test conversions of Pose3Pose3" begin
    dd = convert(PackedPose3Pose3, pp3)
    upd = convert(RoME.Pose3Pose3, dd)


    @test norm(pp3.Zij.μ - upd.Zij.μ) < 1e-10
    @test norm(pp3.Zij.Σ.mat - upd.Zij.Σ.mat) < 1e-8

    packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3}, getData(f2))
    unpackeddata = convert(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose3Pose3}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f2), unpackeddata)
end


odo = SE3(randn(3), convert(Quaternion, so3(0.1*randn(3)) ))
odoc = Pose3Pose3NH( MvNormal(veeEuler(odo),odoCov), [0.5;0.5])

f3 = addFactor!(fg,[:x1;:x2],odoc)


@testset "test conversions of Pose3Pose3NH" begin

    dd = convert(PackedPose3Pose3NH, odoc)
    upd = convert(RoME.Pose3Pose3NH, dd)

    @test norm(odoc.Zij.μ - upd.Zij.μ) < 1e-8
    @test norm(odoc.Zij.Σ.mat - upd.Zij.Σ.mat) < 1e-8
    @test norm(odoc.nullhypothesis.p - upd.nullhypothesis.p) < 1e-8


    packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3NH}, getData(f3))
    unpackeddata = convert(IncrementalInference.FunctionNodeData{GenericWrapParam{RoME.Pose3Pose3NH}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f3), unpackeddata)
end



@testset "test conversions of PartialPriorRollPitchZ" begin

    prpz = PartialPriorRollPitchZ(MvNormal([0.0;0.5],0.1*eye(2)),Normal(3.0,0.5))

    pprpz = convert(PackedPartialPriorRollPitchZ, prpz)
    unp = convert(PartialPriorRollPitchZ, pprpz)

    @test RoME.compare(prpz, unp)

end



@testset "test conversions of PartialPose3XYYaw" begin
    xyy = PartialPose3XYYaw(MvNormal([1.0;2.0;0.5],0.1*eye(3)))

    pxyy = convert(PackedPartialPose3XYYaw, xyy)
    unp = convert(PartialPose3XYYaw, pxyy)

    @test RoME.compare(xyy, unp)
end



@testset "test conversions of PartialPose3XYYawNH" begin

    xyy = PartialPose3XYYawNH(MvNormal([1.0;2.0;0.5],0.1*eye(3)), [0.6;0.4])

    pxyy = convert(PackedPartialPose3XYYawNH, xyy)
    unp = convert(PartialPose3XYYawNH, pxyy)

    @test RoME.compare(xyy, unp)
end


@testset "test PriorPoint2DensityNH" begin

    prpt2 = PriorPoint2DensityNH(kde!(randn(2,100)),[0.25;0.75]  )

    pprpt2 = convert(PackedPriorPoint2DensityNH, prpt2)
    uprpt2 = convert(PriorPoint2DensityNH, pprpt2)

    @test norm(getPoints(prpt2.belief)-getPoints(uprpt2.belief)) < 1e-8

    @test norm(prpt2.nullhypothesis.p-uprpt2.nullhypothesis.p) < 1e-8

end





#
