# test packing functions


using RoME
using Test
using DistributedFactorGraphs
import DistributedFactorGraphs: packVariableNodeData, unpackVariableNodeData


##

@testset "test PriorPoint2" begin
##

global prpt2 = PriorPoint2( MvNormal([0.25;0.75], Matrix(Diagonal([1.0;2.0]))) )

global pprpt2 = convert(PackedPriorPoint2, prpt2)
global uprpt2 = convert(PriorPoint2, pprpt2)

@test norm(prpt2.Z.μ - uprpt2.Z.μ) < 1e-8

@test norm(prpt2.Z.Σ.mat - uprpt2.Z.Σ.mat) < 1e-8

# test backwards compatibility, TODO remove
global prpt2 = PriorPoint2( MvNormal([0.25;0.75], Matrix(Diagonal([1.0;2.0].^2))  ) )

##
end


##

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

##

@testset "test conversions of PriorPose2" begin
global fg
##


dd = convert(PackedPriorPose2, ipp)
upd = convert(RoME.PriorPose2, dd)

@test RoME.compare(ipp, upd)

packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose2}, DFG.getSolverData(f1))
unpackeddata = reconstFactorData(fg, getVariableOrder(f1), IIF.FunctionNodeData{IIF.CommonConvWrapper{RoME.PriorPose2}}, packeddata);
# unpackeddata = convert(IIF.FunctionNodeData{IIF.CommonConvWrapper{RoME.PriorPose2}}, packeddata)

# DFG.getSolverData(f1)
# unpackeddata

@test DFG.compare(DFG.getSolverData(f1), unpackeddata)

# TODO: https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/44
packedv1data = packVariableNodeData(DFG.getSolverData(v1))
upv1data = unpackVariableNodeData(packedv1data)
# packedv1data = convert(IncrementalInference.PackedVariableNodeData, DFG.getSolverData(v1))
# upv1data = convert(IncrementalInference.VariableNodeData, packedv1data)

@test compareAll(DFG.getSolverData(v1), upv1data, skip=[:variableType;:val])
@test all( isapprox.(DFG.getSolverData(v1).val, upv1data.val) )
@test compareFields(DFG.getSolverData(v1).variableType, upv1data.variableType)

##
end


@testset "test conversions of Pose2Pose2" begin
##

global dd = convert(PackedPose2Pose2, ppc)
global upd = convert(RoME.Pose2Pose2, dd)

# TODO -- fix ambiguity in compare function
@test RoME.compare(ppc, upd)

global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Pose2}, DFG.getSolverData(f2))
global unpackeddata = reconstFactorData(fg, getVariableOrder(f2), IIF.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Pose2}}, packeddata)
# global unpackeddata = convert(IIF.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Pose2}}, packeddata)

# TODO -- fix ambibuity in compare function
@test DFG.compare(DFG.getSolverData(f2), unpackeddata)

##
end


@testset "test conversions of Pose2Point2BearingRange" begin
##

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


global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose2Point2BearingRange}, DFG.getSolverData(f3))
global unpackeddata = reconstFactorData(fg, getVariableOrder(f3),IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Point2BearingRange}}, packeddata) # IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose2Point2BearingRange{Normal{Float64},Normal{Float64}}}}

@test ppbr.bearing.μ == unpackeddata.fnc.usrfnc!.bearing.μ
@test ppbr.bearing.σ == unpackeddata.fnc.usrfnc!.bearing.σ

@test ppbr.range.μ == unpackeddata.fnc.usrfnc!.range.μ
@test ppbr.range.σ == unpackeddata.fnc.usrfnc!.range.σ

##
end
##

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

##

@testset "test conversions of PriorPose3" begin

##

global dd = convert(PackedPriorPose3, ipp)
global upd = convert(RoME.PriorPose3, dd)

# @test TransformUtils.compare(ipp.Zi, upd.Zi)
@test norm(ipp.Z.μ - upd.Z.μ) < 1e-10
@test norm(ipp.Z.Σ.mat - upd.Z.Σ.mat) < 1e-8

global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPriorPose3}, DFG.getSolverData(f1))
global unpackeddata = reconstFactorData(fg, getVariableOrder(f1), IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.PriorPose3}}, packeddata)

# TODO -- fix ambibuity in compare function
@test compareAll(DFG.getSolverData(f1), unpackeddata, skip=[:fnc;])
@test_broken compareAll(DFG.getSolverData(f1).fnc, unpackeddata.fnc, skip=[:params;:threadmodel;:cpt;:usrfnc!;:vartypes;:particleidx;:varidx])
@test compareAll(DFG.getSolverData(f1).fnc.usrfnc!, unpackeddata.fnc.usrfnc!, skip=[:Zi;:Z;:p])

# deprecated p
# @test isapprox( DFG.getSolverData(f1).fnc.usrfnc!.p, unpackeddata.fnc.usrfnc!.p)

@test compareAll(DFG.getSolverData(f1).fnc.usrfnc!.Z, unpackeddata.fnc.usrfnc!.Z, skip=[:Σ;])
@test compareAll(DFG.getSolverData(f1).fnc.usrfnc!.Z.Σ, unpackeddata.fnc.usrfnc!.Z.Σ)

@error "not comparing CCW.fnc.cpt, see this test file for details -- pending IIF #825, and DFG #590"
# see https://github.com/JuliaRobotics/DistributedFactorGraphs.jl/issues/590#issuecomment-776838053
# @test compareAll(DFG.getSolverData(f1).fnc.cpt[1], unpackeddata.fnc.cpt[1], skip=[:factormetadata;:activehypo])

# @test compareAll(DFG.getSolverData(f1).fnc.params, unpackeddata.fnc.params)
@warn "threadmodel is obsolete"
# @test compareAll(DFG.getSolverData(f1).fnc.threadmodel, unpackeddata.fnc.threadmodel)

# TODO: Ref above
packedv1data = packVariableNodeData(DFG.getSolverData(v1))
upv1data = unpackVariableNodeData(packedv1data)
# global packedv1data = convert(IncrementalInference.PackedVariableNodeData, DFG.getSolverData(v1))
# global upv1data = convert(IncrementalInference.VariableNodeData, packedv1data)

@test compareAll(DFG.getSolverData(v1), upv1data, skip=[:variableType;:val])
@test all( isapprox.(DFG.getSolverData(v1).val, upv1data.val ) )
@test compareAll(DFG.getSolverData(v1).variableType, upv1data.variableType)

##

end

##

pp3 = Pose3Pose3( MvNormal([25.0;0;0;0;0;0], odoCov) )
v2 = addVariable!(fg, :x2, Pose3)
f2 = addFactor!(fg, [:x1, :x2], pp3)


@testset "test conversions of Pose3Pose3" begin

global dd = convert(PackedPose3Pose3, pp3)
global upd = convert(RoME.Pose3Pose3, dd)


@test norm(pp3.Z.μ - upd.Z.μ) < 1e-10
@test norm(pp3.Z.Σ.mat - upd.Z.Σ.mat) < 1e-8

global packeddata = convert(IncrementalInference.PackedFunctionNodeData{RoME.PackedPose3Pose3}, DFG.getSolverData(f2))
global unpackeddata = reconstFactorData(fg, getVariableOrder(f2), IncrementalInference.FunctionNodeData{IIF.CommonConvWrapper{RoME.Pose3Pose3}}, packeddata)

# TODO -- fix ambibuity in compare function
@test DFG.compare(DFG.getSolverData(f2), unpackeddata)

end


global odo = SE3(randn(3), convert(Quaternion, so3(0.1*randn(3)) ))
global odoc = Pose3Pose3( MvNormal(veeEuler(odo),odoCov))

global f3 = addFactor!(fg,[:x1;:x2],odoc, nullhypo=0.5)



##

@testset "test conversions of PartialPriorRollPitchZ" begin

global prpz = PriorPose3ZRP(MvNormal([0.0;0.5],0.1*diagm([1.0;1])),Normal(3.0,0.5))

global pprpz = convert(PackedPriorPose3ZRP, prpz)
global unp = convert(PriorPose3ZRP, pprpz)

@test RoME.compare(prpz, unp)

end


@testset "test conversions of PartialPose3XYYaw" begin

global xyy = Pose3Pose3XYYaw(
            MvNormal( [1.0;2.0;0.5], 0.1*diagm([1.0;1;1]) ))

global pxyy = convert(PackedPose3Pose3XYYaw, xyy)
global unp = convert(Pose3Pose3XYYaw, pxyy)

@test RoME.compare(xyy, unp)
end

##

@testset "test conversions of Pose3Pose3Rotation" begin

p3rot = Pose3Pose3Rotation(MvNormal( [0.1;0.2;0.3], 0.1*diagm([1.0;1;1]) ))

pac = convert(PackedPose3Pose3Rotation, p3rot)
unp = convert(Pose3Pose3Rotation, p3rot)

@test RoME.compare(p3rot, unp)
end


#
