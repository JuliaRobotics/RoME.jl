# tutorial on conventional 2D SLAM example

# addprocs(3)

# This tutorial shows how to use some of the commonly used factor types
# This tutorial follows from the ContinuousScalar example from IncrementalInference

# @everywhere begin
  using RoME, Distributions
  using IncrementalInference
  using TransformUtils
# end # everywhere


# using Caesar
# backend_config, user_config = standardcloudgraphsetup(addrdict=user_config)
# function initialize!(backend_config,
#                     user_config)
#     println("[Caesar.jl] Setting up factor graph")
#     fg = Caesar.initfg(sessionname=user_config["session"], cloudgraph=backend_config)
#     println("[Caesar.jl] Creating SLAM client/object")
#     return  SLAMWrapper(fg, nothing, 0)
# end


# Start with an empty graph (local dictionary version)
fg = initfg(sessionname="SLAM2D_TUTORIAL")

# also add a PriorPose2 to pin the first pose at a fixed location
addNode!(fg, :x0, Pose2, labels=["POSE"])
addFactor!(fg, [:x0], PriorPose2(zeros(3,1), 0.01*eye(3), [1.0]))
# Prior(MvNormal([0.0;0.0;0], diagm([1.0;1.0;0.01].^2)))

# Drive around in a hexagon
for i in 0:5
  psym = Symbol("x$i")
  nsym = Symbol("x$(i+1)")
  addNode!(fg, nsym, Pose2, labels=["POSE"])
  addFactor!(fg, [psym;nsym], Pose2Pose2(reshape([10.0;0;pi/3],3,1), 0.01*eye(3), [1.0])  )
  # Pose2Pose2_NEW(MvNormal([10.0;0;pi/3], diagm([0.1;0.1;0.1].^2)))
end

# Graphs.plot(fg.g)


ensureAllInitialized!(fg)





tree = wipeBuildNewTree!(fg)
# inferOverTree!(fg, tree)
inferOverTreeR!(fg, tree)





using RoMEPlotting, Gadfly





pl = plotKDE(fg, [:x0; :x1; :x2; :x3; :x4; :x5; :x6]);

Gadfly.draw(PDF("tmpX0123456.pdf", 15cm, 15cm), pl)

@async run(`evince tmpX0123456.pdf`)



# pl = drawPoses(fg)
pl = drawPosesLandms(fg)
Gadfly.draw(PDF("tmpPosesFg.pdf", 15cm, 15cm), pl)
@async run(`evince tmpPosesFg.pdf`)



tree = wipeBuildNewTree!(fg)


@async Graphs.plot(tree.bt)


@time inferOverTree!(fg, tree)


# Graphs.plot(tree.bt)

@async Graphs.plot(fg.g)
