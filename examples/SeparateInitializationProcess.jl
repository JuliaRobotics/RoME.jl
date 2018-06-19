# tutorial on conventional 2D SLAM example

addprocs(3)

# This tutorial shows how to use some of the commonly used factor types
# This tutorial follows from the ContinuousScalar example from IncrementalInference

@everywhere begin
using RoME, Distributions
using IncrementalInference
using TransformUtils

import IncrementalInference: getSample

# struct Pose2 <: IncrementalInference.InferenceVariable
#   dims::Int
#   Pose2() = new(3)
# end
# struct Point2 <: IncrementalInference.InferenceVariable
#   dims::Int
#   Point2() = new(2)
# end


struct Prior{T} <: IncrementalInference.FunctorSingleton where {T <: Distribution}
  z::T
end
getSample(s::Prior, N::Int=1) = (rand(s.z,N), )

# mutable struct Pose2Pose2_NEW{T} <: IncrementalInference.FunctorPairwise where {T <: Distribution}
#   z::T
#   Pose2Pose2_NEW() where {T <: Distribution} = new{T}()
#   Pose2Pose2_NEW(z1::T) where {T <: Distribution} = new{T}(z1)
#   Pose2Pose2_NEW(mean::Vector{Float64}, cov::Array{Float64,2}) where {T <: Distribution} = new{Distributions.MvNormal}(MvNormal(mean, cov))
# end
# getSample(s::Pose2Pose2_NEW, N::Int=1) = (rand(s.z,N), )
# function (s::Pose2Pose2_NEW)(res::Array{Float64},
#             userdata::Union{FactorMetadata,Void},
#             idx::Int,
#             meas::Tuple,
#             wxi::Array{Float64,2},
#             wxj::Array{Float64,2}  )
#   #
#   wXjhat = SE2(wxi[:,idx])*SE2(meas[1][:,idx])
#   jXjhat = SE2(wxj[:,idx]) \ wXjhat
#   se2vee!(res, jXjhat)
#   nothing
# end
# function (s::Pose2Pose2_NEW)(res::Array{Float64},
#       idx::Int,
#       meas::Tuple,
#       wxi::Array{Float64,2},
#       wxj::Array{Float64,2}  )
#   #
#   s(res, nothing, idx, meas, wxi, wxj)
# end
end # everywhere


# Start with an empty graph
fg = initfg(sessionname="SLAM2D_TUTORIAL")


# also add a PriorPose2 to pin the first pose at a fixed location
addNode!(fg, :x0, Pose2, labels=["POSE"])
addFactor!(fg, [:x0], Prior(MvNormal([0.0;0.0;0], diagm([1.0;1.0;0.01].^2))) )

# Drive around in a hexagon
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addNode!(fg, nsym, Pose2)
  addFactor!(fg, [psym;ensureAllInitialized!nsym], Pose2Pose2(MvNormal([10.0;0;pi/3], diagm([0.1;0.1;0.1].^2))))
end

# Graphs.plot(fg.g)

# Add landmarks with Bearing range measurements
addNode!(fg, :l1, Point2, labels=["LANDMARK"])
p2br = Pose2DPoint2DBearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br)


# Add landmarks with Bearing range measurements
p2br2 = Pose2DPoint2DBearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br2)


isInitialized(fg, :l1)


# Replace this call with an initialization service that traverses the entire graph
ensureAllInitialized!(fg)



# get all initialized nodes





using RoMEPlotting, Gadfly



pl = plotKDE(fg, [:x0; :x1; :x2; :x3; :x4; :x5; :x6]);

Gadfly.draw(PDF("tmpX0123456.pdf", 15cm, 15cm), pl)

@async run(`evince tmpX0123456.pdf`)



# pl = drawPoses(fg)
pl = drawPosesLandms(fg)
Gadfly.draw(PDF("tmpPosesFg.pdf", 15cm, 15cm), pl)
@async run(`evince tmpPosesFg.pdf`)



# Perform inference over the entire graph

tree = wipeBuildNewTree!(fg)
# @async Graphs.plot(tree.bt)
@time inferOverTree!(fg, tree)

#
