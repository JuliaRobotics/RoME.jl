
using Test
using RoME


##
@testset "Test covariance weights using Mahalanobis distance (for parametric solve)" begin
## sanity check

RES = Float64[]

SW = -0.1:0.0001:0.0
for x_test in SW
  # stronger
  x1 = +0.1
  c1 =  0.01
  # weaker
  x2 = -0.1
  c2 =  0.03
  # Mahalanobis distance
  dx = x_test .- [x1;x2]
  iS = [c1 0; 0 c2]
  res = dx' * iS * dx
  push!(RES, res)
end

@test isapprox( -0.05, SW[ findfirst(s-> s == minimum(RES), RES) ]; atol=1e-3)

##
end


@testset "Test parametric solve covariances, Pose2" begin
## Test Pose2

fg = initfg()

addVariable!(fg, :x0, Pose2)
addVariable!(fg, :x1, Pose2)

addFactor!(fg, [:x0], PriorPose2(MvNormal([0.,0,0], diagm([0.1,0.1,0.01].^2))))

addFactor!(fg, [:x0,:x1], Pose2Pose2(MvNormal([1.1,0,0], diagm([0.1,0.1,0.01].^2))))
addFactor!(fg, [:x0,:x1], Pose2Pose2(MvNormal([0.9,0,0], diagm([sqrt(0.03),0.1,0.01].^2))))

##

IIF.autoinitParametric!(fg)
# IIF.solveGraph!(fg)

@test isapprox( [0;0;0.], getPPESuggested(fg, :x0, :parametric); atol=1e-4 )
@test isapprox( [1.05;0;0], getPPESuggested(fg, :x1, :parametric); atol=1e-4 )

##
end
##

# using GLMakie
# lines(SW, RES)

##