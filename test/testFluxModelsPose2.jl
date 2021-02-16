# test FluxModelsPose2Pose2

using Test
using LinearAlgebra
using Flux
using RoME

@testset "Basic funtionality tests on FluxModelsPose2Pose2" begin


##

mdls = [RoME.buildPose2OdoNN_01_FromElements( randn(4,8),
                                              randn(8),
                                              randn(8,48),
                                              randn(8),
                                              randn(2,8),
                                              randn(2)) for i in 1:10];


##

# start with a basic factor graph

mvnNaive = MvNormal(zeros(3), diagm([1.0;1.0;0.01]))

jvd = zeros(25,4)
pp = FluxModelsPose2Pose2(mdls, jvd, mvnNaive, 0.5)

##


fg = generateCanonicalFG_ZeroPose2()
addVariable!(fg, :x1, Pose2)

addFactor!(fg, [:x0;:x1], pp)

##

pts = approxConv(fg, :x0x1f1, :x1)

#


end



