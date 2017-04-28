
using Caesar
using RoME, Distributions
import IncrementalInference: getSample

using TransformUtils


type HipJoint <: FunctorPairwise
  Zij::Normal
end
function getSample(el::HipJoint, N=1)
  return (rand(el.Zij, N),)
end
function (el::HipJoint)(res, idx, meas, xi, xj)
  base = SE3(xi[1:3,idx],Euler(xi[4:6,idx]))
  h=meas[1][idx]
  hip = SE3([0,0,0.0], Euler(0,0,h))
  torso = ( base ⊕ hip ) ⊕ SE3([0,0,0.5],SO3(0))

  delta = torso\SE3(xj[1:3,idx],Euler(xj[4:6,idx]))

  res[1:6] = veeEuler(delta)
end




N=100
fg = initfg()

# base
addNode!(fg, :x1, dims=6)
pos = PriorPose3(MvNormal(zeros(6),0.000001*eye(6)))
addFactor!(fg, [:x1], pos) # base
initializeNode!(fg, :x1)

# torso
addNode!(fg, :x2, dims=6)
hip = HipJoint(Normal(0,0.1))
addFactor!(fg, [:x1, :x2], hip) # hio
initializeNode!(fg, :x2)




using IncrementalInference

pts = evalFactor2(fg, fg.g.vertices[fg.fIDs[:x1x2]], fg.IDs[:x2])

@show Base.mean(pts,2)



# setup visualization process and default drawings
vis = startdefaultvisualization()


visualizeallposes!(vis, fg, drawtype=:max)



solveandvisualize(fg, vis)

Graphs.plot(fg.g)

visualizeDensityMesh!(vis, fg, :x2)


plotPose3Pairs(fg, :x2)




















type ShoulderJoint <: FunctorPairwise
  Zij::Normal
end
function getSample(el::ShoulderJoint, N=1)
  return (rand(el.Zij, N),)
end
function (el::ShoulderJoint)(res, idx, meas, xi, xj)
  hip = SE3(xi[1:3,idx],Euler(xi[4:6,idx]))
  s=el.Zij.μ
  hip = SE3([0,0,0.0], Euler(s,0,0))
  arm = ( xbase ⊕ hip ) ⊕ SE3([0,0,0.4],SO3(0))

  delta = arm\SE3(xj[1:3,idx],Euler(xj[4:6,idx]))

  res[1:6] = veeEuler(delta)
end
