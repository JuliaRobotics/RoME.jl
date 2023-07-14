# using Revise
using RoME
using Manifolds
using StaticArrays
using Test

##
@testset "test SE(3) coordinates to homography and back" begin
##

M = getManifold(Pose3)

@test M == Manifolds.SpecialEuclidean(3)

C = 0.2*randn(6)
H = coordinates_to_homography(M, C)
C_ = homography_to_coordinates(M, H)

@test isapprox(C, C_)

##
end

@testset "Test (partial) priors on Pose3" begin
## there was a problem with partials and staticarrays where asigning PriorPoint3 to Pose3 stopped working, here is an explicit test

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :x0, Pose3)

mv3 = MvNormal(SA[0.0; 0.0; 0.0], SMatrix{3,3}(1,0,0,0,1,0,0,0,1.))
f = addFactor!(fg, [:x0;], PriorPoint3(mv3))

approxConvBelief(fg, getLabel(f), :x0)

##
end

##
@testset "Test Basic Pose3 :parametric and :default" begin
fg = initfg()

prior_distribution = PriorPose3( MvNormal([0., 0, 0, 0, -pi/4, 0], diagm([0.1, 0.1, 0.1, 0.01, 0.01, 0.01]).^2))

addVariable!(fg, :x0, Pose3)
prior = addFactor!(fg, [:x0], prior_distribution) 

odo_distribution = MvNormal([sqrt(2), 0, 0, 0, 0, pi/2], diagm([0.1, 0.1, 0.1, 0.01, 0.01, 0.01].^2))

for i=1:4
    f = Symbol("x",i-1)
    t = Symbol("x",i)
    addVariable!(fg, t, Pose3)
    addFactor!(fg, [f, t], Pose3Pose3(odo_distribution))
end

initAll!(fg)

@time r = IIF.solveGraphParametric!(fg)

M = getManifold(Pose3)

p0 = getVal(fg, :x0, solveKey=:parametric)[1]
p1 = getVal(fg, :x1, solveKey=:parametric)[1]
p2 = getVal(fg, :x2, solveKey=:parametric)[1]
p3 = getVal(fg, :x3, solveKey=:parametric)[1]
p4 = getVal(fg, :x4, solveKey=:parametric)[1]

@test isapprox(M, p0, p4, atol=0.001)

np0 = mean(M, getVal(fg, :x0))
np1 = mean(M, getVal(fg, :x1))
np2 = mean(M, getVal(fg, :x2))
np3 = mean(M, getVal(fg, :x3))
np4 = mean(M, getVal(fg, :x4))

@test isapprox(M, p0, np0, atol=0.1)
@test isapprox(M, p1, np1, atol=0.1)
@test isapprox(M, p2, np2, atol=0.1)
@test isapprox(M, p3, np3, atol=0.1)
@test isapprox(M, p4, np4, atol=0.1)

end
