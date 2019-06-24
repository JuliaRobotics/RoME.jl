# test packing functions


# using Distributions
# using KernelDensityEstimate
# using TransformUtils
using RoME
using Test
using DistributedFactorGraphs
import DistributedFactorGraphs: pack, unpack



@testset "test PriorPoint2" begin

  global prpt2 = PriorPoint2( MvNormal([0.25;0.75], Matrix(Diagonal([1.0;2.0]))) )

  global pprpt2 = convert(PackedPriorPoint2, prpt2)
  global uprpt2 = convert(PriorPoint2, pprpt2)

  @test norm(prpt2.Z.μ - uprpt2.Z.μ) < 1e-8

  @test norm(prpt2.Z.Σ.mat - uprpt2.Z.Σ.mat) < 1e-8

  # test backwards compatibility, TODO remove
  global prpt2 = PriorPoint2( MvNormal([0.25;0.75], Matrix(Diagonal([1.0;2.0].^2))  ) )

end




N = 100
fg = initfg()


initCov = Matrix(Diagonal([0.03;0.03;0.001]))
odoCov = Matrix(Diagonal([3.0;3.0;0.01]))

# Some starting position
v1 = addVariable!(fg, :x1, Pose2, N=N)
ipp = PriorPose2(MvNormal(zeros(3), initCov))
f1  = addFactor!(fg,[v1], ipp)

# and a second pose
v2 = addVariable!(fg, :x2, Pose2, N=N)
ppc = Pose2Pose2( MvNormal([50.0;0.0;pi/2], odoCov) )
f2 = addFactor!(fg, [:x1;:x2], ppc)



@testset "test conversions of PriorPose2" begin
global fg

dd = convert(PackedPriorPose2, ipp)
upd = convert(RoME.PriorPose2, dd)

@test RoME.compare(ipp, upd) # temp use of RoME.compare

packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose2}, getData(f1))
unpackeddata = convert(IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.PriorPose2}}, packeddata)

# getData(f1)
# unpackeddata

@test RoME.compare(getData(f1), unpackeddata) # temp use of RoME.compare
# TODO -- what what??
# setSerializationNamespace!("Main" => Main)
# TODO: https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/44
packedv1data = pack(fg, getData(v1))
upv1data = unpack(fg, packedv1data)
# packedv1data = convert(IncrementalInference.PackedVariableNodeData, getData(v1))
# upv1data = convert(IncrementalInference.VariableNodeData, packedv1data)

@test compareAll(getData(v1), upv1data, skip=[:softtype;])
@test compareFields(getData(v1).softtype, upv1data.softtype)

end

@testset "test conversions of Pose2Pose2" begin

    global dd = convert(PackedPose2Pose2, ppc)
    global upd = convert(RoME.Pose2Pose2, dd)

    # TODO -- fix ambiguity in compare function
    @test RoME.compare(ppc, upd)

    global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Pose2}, getData(f2))
    global unpackeddata = convert(IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Pose2}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f2), unpackeddata)
end


@testset "test conversions of Pose2Point2BearingRange" begin
    # and a second pose
    global v3 = addVariable!(fg, :l1, Point2, N=N)
    # ppc = Pose2Point2BearingRange([50.0;0.0;pi/2], 0.01*Matrix{Float64}(LinearAlgebra.I, 2,2), [1.0])
    global ppbr = Pose2Point2BearingRange(
                  Normal(0.0, 0.005 ),
                  Normal(50, 0.5) )
    global f3 = addFactor!(fg, [:x2;:l1], ppbr)

    global dd = convert(PackedPose2Point2BearingRange, ppbr)
    global upd = convert(
            RoME.Pose2Point2BearingRange,
            dd
          )


    global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Point2BearingRange}, getData(f3))
    global unpackeddata = convert(IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Point2BearingRange}}, packeddata) # IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Point2BearingRange{Normal{Float64},Normal{Float64}}}}

    @test ppbr.bearing.μ == unpackeddata.fnc.usrfnc!.bearing.μ
    @test ppbr.bearing.σ == unpackeddata.fnc.usrfnc!.bearing.σ

    @test ppbr.range.μ == unpackeddata.fnc.usrfnc!.range.μ
    @test ppbr.range.σ == unpackeddata.fnc.usrfnc!.range.σ
end

# parameters
global N = 300
global initCov = Matrix{Float64}(LinearAlgebra.I, 6,6)
[initCov[i,i] = 0.01 for i in 4:6];
global odoCov = deepcopy(initCov)

# start new factor graph
global fg = initfg()

global v1 = addVariable!(fg, :x1, Pose3, N=N) #  0.1*randn(6,N)
global ipp = PriorPose3(MvNormal(zeros(6), initCov) )
global f1  = addFactor!(fg,[v1], ipp)


@testset "test conversions of PriorPose3" begin

    global dd = convert(PackedPriorPose3, ipp)
    global upd = convert(RoME.PriorPose3, dd)

    # @test TransformUtils.compare(ipp.Zi, upd.Zi)
    @test norm(ipp.Zi.μ - upd.Zi.μ) < 1e-10
    @test norm(ipp.Zi.Σ.mat - upd.Zi.Σ.mat) < 1e-8

    global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose3}, getData(f1))
    global unpackeddata = convert(IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.PriorPose3}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test compareAll(getData(f1), unpackeddata, skip=[:fnc;])
    @test compareAll(getData(f1).fnc, unpackeddata.fnc, skip=[:params;:threadmodel;:cpt;:usrfnc!])
    @test compareAll(getData(f1).fnc.usrfnc!, unpackeddata.fnc.usrfnc!, skip=[:Zi;])
    @test compareAll(getData(f1).fnc.usrfnc!.Zi, unpackeddata.fnc.usrfnc!.Zi, skip=[:Σ;])
    @test compareAll(getData(f1).fnc.usrfnc!.Zi.Σ, unpackeddata.fnc.usrfnc!.Zi.Σ)
    @test compareAll(getData(f1).fnc.cpt[1], unpackeddata.fnc.cpt[1], skip=[:factormetadata;:activehypo])
    # @test compareAll(getData(f1).fnc.params, unpackeddata.fnc.params)
    @warn "threadmodel is not defined, fix with DFG"
    # @test compareAll(getData(f1).fnc.threadmodel, unpackeddata.fnc.threadmodel)

    # TODO: Ref above
    packedv1data = pack(fg, getData(v1))
    upv1data = unpack(fg, packedv1data)
    # global packedv1data = convert(IncrementalInference.PackedVariableNodeData, getData(v1))
    # global upv1data = convert(IncrementalInference.VariableNodeData, packedv1data)

    @test compareAll(getData(v1), upv1data, skip=[:softtype;])
    @test compareAll(getData(v1).softtype, upv1data.softtype)

end



pp3 = Pose3Pose3( MvNormal([25.0;0;0;0;0;0], odoCov) )
v2 = addVariable!(fg, :x2, Pose3)
f2 = addFactor!(fg, [:x1, :x2], pp3)


@testset "test conversions of Pose3Pose3" begin
    global dd = convert(PackedPose3Pose3, pp3)
    global upd = convert(RoME.Pose3Pose3, dd)


    @test norm(pp3.Zij.μ - upd.Zij.μ) < 1e-10
    @test norm(pp3.Zij.Σ.mat - upd.Zij.Σ.mat) < 1e-8

    global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3}, getData(f2))
    global unpackeddata = convert(IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose3Pose3}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f2), unpackeddata)
end


global odo = SE3(randn(3), convert(Quaternion, so3(0.1*randn(3)) ))
global odoc = Pose3Pose3NH( MvNormal(veeEuler(odo),odoCov), [0.5;0.5])

global f3 = addFactor!(fg,[:x1;:x2],odoc)


@testset "test conversions of Pose3Pose3NH" begin

    global dd = convert(PackedPose3Pose3NH, odoc)
    global upd = convert(RoME.Pose3Pose3NH, dd)

    @test norm(odoc.Zij.μ - upd.Zij.μ) < 1e-8
    @test norm(odoc.Zij.Σ.mat - upd.Zij.Σ.mat) < 1e-8
    @test norm(odoc.nullhypothesis.p - upd.nullhypothesis.p) < 1e-8


    global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3NH}, getData(f3))
    global unpackeddata = convert(IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose3Pose3NH}}, packeddata)

    # TODO -- fix ambibuity in compare function
    @test IncrementalInference.compare(getData(f3), unpackeddata)
end



@testset "test conversions of PartialPriorRollPitchZ" begin

    global prpz = PartialPriorRollPitchZ(MvNormal([0.0;0.5],0.1*Matrix{Float64}(LinearAlgebra.I, 2,2)),Normal(3.0,0.5))

    global pprpz = convert(PackedPartialPriorRollPitchZ, prpz)
    global unp = convert(PartialPriorRollPitchZ, pprpz)

    @test RoME.compare(prpz, unp)

end


@testset "test conversions of PartialPose3XYYaw" begin
    global xyy = PartialPose3XYYaw(
             MvNormal( [1.0;2.0],
                        0.1*Matrix{Float64}(LinearAlgebra.I, 2,2) ),
             Normal(0.5, 0.1)
           )

    global pxyy = convert(PackedPartialPose3XYYaw, xyy)
    global unp = convert(PartialPose3XYYaw, pxyy)

    @test RoME.compare(xyy, unp)
end



@testset "test conversions of PartialPose3XYYawNH" begin

    global xyy = PartialPose3XYYawNH(MvNormal([1.0;2.0;0.5],0.1*Matrix{Float64}(LinearAlgebra.I, 3,3)), [0.6;0.4])

    global pxyy = convert(PackedPartialPose3XYYawNH, xyy)
    global unp = convert(PartialPose3XYYawNH, pxyy)

    @test RoME.compare(xyy, unp)
end


@testset "test PriorPoint2DensityNH" begin

    global prpt2 = PriorPoint2DensityNH(kde!(randn(2,100)),[0.25;0.75]  )

    global pprpt2 = convert(PackedPriorPoint2DensityNH, prpt2)
    global uprpt2 = convert(PriorPoint2DensityNH, pprpt2)

    @test norm(getPoints(prpt2.belief)-getPoints(uprpt2.belief)) < 1e-8

    @test norm(prpt2.nullhypothesis.p-uprpt2.nullhypothesis.p) < 1e-8

end





#
