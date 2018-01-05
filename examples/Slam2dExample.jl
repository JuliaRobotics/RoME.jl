# tutorial on conventional 2D SLAM example

# This tutorial shows how to use some of the commonly used factor types
# This tutorial follows from the ContinuousScalar example from IncrementalInference

using RoME, Distributions
using IncrementalInference
using TransformUtils

import IncrementalInference: getSample

struct Pose2 <: IncrementalInference.InferenceVariable
  dims::Int
  Pose2() = new(3)
end

struct Prior{T} <: IncrementalInference.FunctorSingleton where {T <: Distribution}
  z::T
end
getSample(s::Prior, N::Int=1) = (rand(s.z,N), )

struct Pose2Pose2_NEW{T} <: IncrementalInference.FunctorPairwise where {T <: Distribution}
  z::T
end
getSample(s::Pose2Pose2_NEW, N::Int=1) = (rand(s.z,N), )
function (s::Pose2Pose2_NEW)(res::Array{Float64},
      idx::Int,
      meas::Tuple,
      wxi::Array{Float64,2},
      wxj::Array{Float64,2}  )
  # res[1] = meas[1][idx] - (X2[1,idx] - X1[1,idx])
  wXjhat = SE2(wxi[:,idx])*SE2(meas[1][:,idx]) #*SE2(pp2.Zij[:,1])*SE2(meas[1][:,idx])
  jXjhat = SE2(wxj[:,idx]) \ wXjhat
  se2vee!(res, jXjhat)
  nothing
end

# Start with an empty graph
fg = initfg(sessionname="SLAM2D_TUTORIAL")


addNode!(fg, :x0, Pose2, labels=["POSE"])

# also add a PriorPose2 to pin the first pose at a fixed location
addFactor!(fg, [:x0], Prior(MvNormal(zeros(3), diagm([1.0;1.0;0.1]))) )


addNode!(fg, :x1, Pose2, labels=["POSE"])

addFactor!(fg, [:x0;:x1], Pose2Pose2_NEW(MvNormal([10.0;0;pi/2], diagm([0.1;0.1;0.03].^2)))
)

ls(fg, :x1)

getSample(getData(getVert(fg, :x0x1f1, nt=:fnc)).fnc.usrfnc!,5)[1]

# Graphs.plot(fg.g)

using RoMEPlotting, Gadfly

pl = plotKDE(fg, :x0);


Gadfly.draw(PDF("tmpX0.pdf", 15cm, 15cm), pl)

run(`evince tmpX0.pdf`)


Graphs.plot(fg.g)


ensureAllInitialized!(fg)



pl = plotKDE(fg, [:x0; :x1]);

Gadfly.draw(PDF("tmpX01.pdf", 15cm, 15cm), pl)


#
