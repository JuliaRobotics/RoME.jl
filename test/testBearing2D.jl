# test bearing

using RoME
using Test
using TensorCast
using DistributedFactorGraphs
using Manifolds: hat

##

@testset "Testing Bearing2D factor" begin
##
M = SpecialEuclidean(2)
ϵ = identity_element(M)
ps = [exp(M, ϵ,  hat(M, ϵ, [0.,0,0]))]
push!(ps, exp(M, ϵ,  hat(M, ϵ, [5.,0,0])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [10.,0,0])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [10.,5,0])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [10.,10,0])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [5.,10,0])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [0.,10,0])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [0.,5,0])))

push!(ps, exp(M, ϵ,  hat(M, ϵ, [0.,0,pi/4])))
push!(ps, exp(M, ϵ,  hat(M, ϵ, [0.,0,-pi/4])))

push!(ps, exp(M, ϵ,  hat(M, ϵ, [1.,2,0]))) # [4,3,0]
# push!(ps, exp(M, ϵ,  hat(M, ϵ, [1.,2,pi/2]))) # [3,-4,-pi/2]
# push!(ps, exp(M, ϵ,  hat(M, ϵ, [1.,2,pi]))) # [-4,-3,-pi]
# push!(ps, exp(M, ϵ,  hat(M, ϵ, [1.,2,-pi]))) # [-3,4,pi]

rs = [0, -pi/4, -pi/2, -3pi/4, -pi, 3pi/4, pi/2, pi/4, pi/4, -pi/4]

push!(rs, pi/4 - atan(3,4))

q = [5., 5]
m = [pi/4]

f = Pose2Point2Bearing(Normal(pi/4,0.05))

for (i,p) in enumerate(ps)
    # @show p
    r = calcFactorResidualTemporary(f, (Pose2, Point2), m, (p, q))
    # @warn "comp: $(isapprox(r, rs[i], atol=0.001))" a=r b=rs[i]
    @test isapprox(r, rs[i], atol=0.001) 
end

@warn "Bearing2D, must still test factor gradients, which will also verify the sign on residual function calculations"
@test_broken false

# test sign
f = Pose2Point2Bearing(Normal(pi/2,0.001))
xi = ProductRepr([0.,0], [1. 0; 0 1])
xj = [1.,1]
res = calcFactorResidualTemporary(f, (Pose2, Point2), [], (xi, xj))
@test isapprox(res, pi/4, atol=0.1) # FIXME, confirm the sign

# test -pi +pi case
f = Pose2Point2Bearing(Normal(pi,0.001))
xi = ProductRepr([0.,0], [1. 0; 0 1])
xj = [-1, -0.001]
res = calcFactorResidualTemporary(f, (Pose2, Point2), [pi], (xi, xj))
@test isapprox(res, -0.001, atol=1e-3)

##
end


@testset "Simple Bearing2D test to give narrow vs broad belief posteriors" begin
## new factor graph

fg = initfg()

addVariable!(fg, :x1, Pose2)
addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x1], PriorPose2(MvNormal([10.0; 0; 0], diagm([0.01;0.01;0.001].^2))))
addFactor!(fg, [:x2], PriorPose2(MvNormal([0.0;10.0; 0], diagm([0.01;0.01;0.001].^2))))

# limit the scope nearby origin, but test was written because it was found with Bearing
# that many points would suddenly land very far away from the origin, even with the prior.
addVariable!(fg, :l1, Point2)
addFactor!(fg, [:l1], PriorPoint2(MvNormal([0.,0], diagm([10.,10].^2))))

initAll!(fg)

addFactor!(fg, [:x1;:l1], Pose2Point2Bearing(Normal(pi, 0.001)))
addFactor!(fg, [:x2;:l1], Pose2Point2Bearing(Normal(-pi/2,0.001)))

points = approxConv(fg, :x1l1f1, :l1)
@cast pts[j,i] := points[i][j]
#x1l1 "partial" constrian on x-axis 
#FIXME very pessimistic tests, should be way better
@test 30 < sum(abs.(pts[2,:]) .< 100)

L1 = approxConvBelief(fg, :x2l1f1, :l1)
points = getPoints(L1, false)
@cast pts[j,i] := points[i][j]
@test 30 < sum(abs.(pts[1,:]) .< 100)

##

solveTree!(fg)

points = getPoints(getBelief(fg, :l1))
@cast pts[j,i] := points[i][j]
pts = collect(pts)
#FIXME check test after Bearing2D is fixed
@test_broken all([80,80] .< sum(abs.(pts) .< [10,10],dims=2) )

##
end

@testset "Triangulation test in 2D, 3 beacons" begin
##

# noise models
lmp_noise = diagm([0.01;0.01].^2)

# new factor graph
fg = initfg()

# landmarks
addVariable!(fg, :l1, Point2)
addVariable!(fg, :l2, Point2)
addVariable!(fg, :l3, Point2)

addFactor!(fg, [:l1], PriorPoint2(MvNormal([10.0;1.0],lmp_noise)), graphinit=false)
addFactor!(fg, [:l2], PriorPoint2(MvNormal([10.0+sqrt(3)/2;-0.5],lmp_noise)), graphinit=false)
addFactor!(fg, [:l3], PriorPoint2(MvNormal([10.0-sqrt(3)/2;-0.5],lmp_noise)), graphinit=false)

# pose
addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x1;:l1], Pose2Point2Bearing(Normal(pi/2,0.05)), graphinit=false)
addFactor!(fg, [:x1;:l2], Pose2Point2Bearing(Normal(-pi/6,0.05)), graphinit=false)
addFactor!(fg, [:x1;:l3], Pose2Point2Bearing(Normal(-pi+pi/6,0.05)), graphinit=false)

# #TODO remove easier just for tests
# # landmarks
# addFactor!(fg, [:l1], PriorPoint2(MvNormal([10.0;0.0],lmp_noise)))
# addFactor!(fg, [:l2], PriorPoint2(MvNormal([0.0;10.0],lmp_noise)))
# addFactor!(fg, [:l3], PriorPoint2(MvNormal([10.0,10],lmp_noise)))
# # pose
# addVariable!(fg, :x1, Pose2)
# addFactor!(fg, [:x1;:l1], Pose2Point2Bearing(Normal(0.0,0.05)), graphinit=false)
# addFactor!(fg, [:x1;:l2], Pose2Point2Bearing(Normal(pi/2,0.05)), graphinit=false)
# addFactor!(fg, [:x1;:l3], Pose2Point2Bearing(Normal(pi/4,0.05)), graphinit=false)


# initVariable!(fg, :x1, [30.0*randn(2,100);randn(1,100)])


## solve

# getSolverParams(fg).drawtree = true
# getSolverParams(fg).showtree = true


tree = solveTree!(fg)


## Look at results
#
# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
# plotSLAM2D(fg, spscale=0.2) # |> PNG("/home/dehann/Downloads/triangulation.png",20cm,15cm)
# plotPose(fg, :x1, scale=0.1) # |> PNG("/home/dehann/Downloads/triangulationPose.png",20cm,15cm)
# plotFactor(fg,  lsf(fg, Pose2Point2Bearing)[1])
# plotLocalProduct(fg, :x1)

## complete the unit test

points = getPoints(getBelief(fg, :x1))
mean(SpecialEuclidean(2), points)

@cast pts[j,i] := getCoordinates.(Pose2, points)[i][j]
pts = collect(pts)
pts[1,:] .-= 10.0

N = size(pts,2)

@test_broken 0.8*N < sum(sqrt.(sum(pts[1:2,:].^2,dims=1)) .< 0.3)

@test_broken 0.8*N < sum(abs.(pts[3,:]) .< 0.1)

#

##
end




@testset "Triangulation test in 2D (opposite), 3 beacons" begin
##

# noise models
lmp_noise = Matrix(Diagonal([0.01;0.01].^2))

# new factor graph
fg = initfg()

# landmarks
addVariable!(fg, :l1, Point2)
addVariable!(fg, :l2, Point2)
addVariable!(fg, :l3, Point2)

addFactor!(fg, [:l1], PriorPoint2(MvNormal([-10.0;1.0-10.0],lmp_noise)), graphinit=false)
addFactor!(fg, [:l2], PriorPoint2(MvNormal([-10.0+sqrt(3)/2;-0.5-10.0],lmp_noise)), graphinit=false)
addFactor!(fg, [:l3], PriorPoint2(MvNormal([-10.0-sqrt(3)/2;-0.5-10.0],lmp_noise)), graphinit=false)

# pose
addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x1;:l1], Pose2Point2Bearing(Normal(pi/2,0.05)), graphinit=false)
addFactor!(fg, [:x1;:l2], Pose2Point2Bearing(Normal(-pi/6,0.05)), graphinit=false)
addFactor!(fg, [:x1;:l3], Pose2Point2Bearing(Normal(-pi+pi/6,0.05)), graphinit=false)


# initVariable!(fg, :x1, [0.01*randn(2,100);-randn(1,100)])
# initVariable!(fg, :x1, [30.0*randn(2,100);randn(1,100)])

# drawGraph(fg)

tree = solveTree!(fg)

points = getPoints(getBelief(fg, :x1))
mean(SpecialEuclidean(2), points)

@cast pts[j,i] := getCoordinates.(Pose2, points)[i][j]
pts = collect(pts)

pts[1,:] .+= 10.0
pts[2,:] .+= 10.0

N = size(pts,2)

@test_broken 0.8*N < sum(sqrt.(sum(pts[1:2,:].^2,dims=1)) .< 0.3)

@test_broken 0.8*N < sum(abs.(pts[3,:]) .< 0.1)

#
# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)
# drawPosesLandms(fg, spscale=0.2) # |> PNG("/home/dehann/Downloads/triangulation.png",20cm,15cm)
# plotPose(fg, :x1, scale=0.1) # |> PNG("/home/dehann/Downloads/triangulationPose.png",20cm,15cm)
#
# plotLocalProduct(fg, :x1)
#
# plotFactor(fg,  lsf(fg, Pose2Point2Bearing)[1])

##
end





#
