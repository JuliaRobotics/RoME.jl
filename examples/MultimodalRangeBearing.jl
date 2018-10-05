using RoME, Distributions
# using RoMEPlotting

import IncrementalInference: getSample

mutable struct NorthSouthPartial{T} <: FunctorSingleton
  Z::T
  partial::Tuple
  NorthSouthPartial() = new()
  NorthSouthPartial(Z::D) where {D <: Distribution} = new{D}(Z, (2,))
end
getSample(ns::NorthSouthPartial, N=1) = (reshape(rand(ns.Z, N),1,N),)


# Start with an empty graph
fg = initfg(sessionname="MULTIMODAL_2D_TUTORIAL")


# Add landmarks with Bearing range measurements
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
addNode!(fg, :l2, Point2, labels=["LANDMARK"])

addFactor!(fg, [:l1], Prior(MvNormal([10.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) )
addFactor!(fg, [:l2], Prior(MvNormal([30.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))) )

setVal!(getVert(fg, :l2), zeros(2,100))

addNode!(fg, :x0, Pose2)
# addFactor!(fg, [:x0], Prior(MvNormal([0.0;0.0;0], Matrix(Diagonal([1.0;1.0;0.01].^2)))) )


p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1; :l2], p2br, multihypo=[1.0; 0.5; 0.5])

addFactor!(fg, [:x0;], NorthSouthPartial(Normal(0,1.0)))

# writeGraphPdf(fg)

ensureAllInitialized!(fg)

tree = wipeBuildNewTree!(fg)

[inferOverTreeR!(fg, tree, N=100) for i in 1:4]




# L1 = getVal(fg, :l1)
# L2 = getVal(fg, :l2)
X0 = getVertKDE(fg, :x0)

# plot(X0, dims=[1;2])


X0pts = getPoints(X0)


@test sum([0 < sum(-20 .< X0pts[1,:] .< 0);
           0 < sum(  0 .< X0pts[1,:] .< 20);
           0 < sum( 20 .< X0pts[1,:] .< 40);
           0 < sum( 40 .< X0pts[1,:] .< 60)]) > 2



# # Drive around in a hexagon
# for i in 0:5
#   psym = Symbol("x$i")
#   nsym = Symbol("x$(i+1)")
#   addNode!(fg, nsym, Pose2, labels=["POSE"])
#   addFactor!(fg, [psym;nsym], Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2)))))
# end
