using RoME
using Test
# using RoMEPlotting, Distributions

import IncrementalInference: getSample

mutable struct NorthSouthPartial{T} <: AbstractPrior
  Z::T
  partial::Tuple{Int}
end

NorthSouthPartial(Z::D) where {D <: IIF.SamplableBelief} = NorthSouthPartial(Z, (2,))

getSample(cfo::CalcFactor{<:NorthSouthPartial}, N::Int=1) = ([[rand(cfo.factor.Z)] for _=1:N],)

##

@testset "test standard multimodal range bearing factor setup..." begin

##

# Start with an empty graph
global N = 100
global fg = initfg(sessionname="MULTIMODAL_2D_TUTORIAL")
# fg.solverParams.attemptGradients=false
# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
addVariable!(fg, :l2, Point2, tags=[:LANDMARK;])

addFactor!(fg, [:l1], PriorPoint2(MvNormal([10.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) )
addFactor!(fg, [:l2], PriorPoint2(MvNormal([30.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) )

val, = predictbelief(fg, :l1, [:l1f1;])
setVal!(fg, :l1, val)

val, = predictbelief(fg, :l2, [:l2f1;])
setVal!(fg, :l2, val)


addVariable!(fg, :x0, Pose2)
# addFactor!(fg, [:x0], Prior(MvNormal([0.0;0.0;0], Matrix(Diagonal([1.0;1.0;0.01].^2)))) )
addFactor!(fg, [:x0;], NorthSouthPartial(Normal(0,1.0)))


global p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1; :l2], p2br, multihypo=[1.0; 0.5; 0.5])


predictbelief(fg, :x0, :) #ls(fg, :x0))

##

end


@testset "test multimodal bearing range factors calculate pose position properly..." begin

##

solveTree!(fg);


_X0pts = getBelief(fg, :x0) |> getPoints
@cast cX0pts[j,i] := getCoordinates.(Pose2, _X0pts)[i][j]
global X0pts = copy(cX0pts)

@test size(X0pts,1) == 3
@test size(X0pts,2) == N

# find at least two of the four possible modes
@test sum([0 < sum(-20 .< X0pts[1,:] .< 0);
           0 < sum(  0 .< X0pts[1,:] .< 20);
           0 < sum( 20 .< X0pts[1,:] .< 40);
           0 < sum( 40 .< X0pts[1,:] .< 60)]) >= 2

##

end


@testset "test multimodal landmark locations are computed correclty..." begin

##

# Start with an empty graph
fg = initfg(sessionname="MULTIMODAL_2D_TUTORIAL")
# getSolverParams(fg).spreadNH = 10.0

addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal([0.0;0.0;0], diagm([1.0;1.0;0.01].^2))) )

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
addFactor!(fg, [:l1], PriorPoint2(MvNormal([40.0;0.0], diagm([1.0;1.0].^2))) )

addVariable!(fg, :l2, Point2, tags=[:LANDMARK;])
addFactor!(fg, [:l2;], NorthSouthPartial(Normal(0,1.0)))
# addFactor!(fg, [:l2], PriorPose2(MvNormal([30.0;0.0], diagm([1.0;1.0].^2))) )

p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1; :l2], p2br, multihypo=[1.0; 0.5; 0.5])

# solve the graph
solveTree!(fg);

# check 


# check modes on L2
L2 = getBelief(fg, :l2)

@cast L2pts[j,i] := getPoints(L2)[i][j]

@test size(L2pts,1) == 2
@test size(L2pts,2) == N


# some likelihood that L2 is around +20
mask = 10 .< L2pts[1,:] .< 30
numM1 = sum(mask)
@test 20 < numM1 < 70

# should also have likelihood of being elsewhere
imask = xor.(mask, 1)
numM2 = sum(imask)
@test 20 < numM2 < 70

##

# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm, 20cm)

##

# plotKDE(fg, :l1, dims=[1])
# plotKDE(L2, dims=[1])
# plotSLAM2D(fg)

# plotLocalProduct(fg, :l1, levels=10)
# plotLocalProduct(fg, :l2, levels=10)

##

end


