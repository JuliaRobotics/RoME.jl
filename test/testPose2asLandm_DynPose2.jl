# build test case for using Pose2Pose2 to landmarks -- to help AprilTag

using RoME  #, Distributions
using RoMEPlotting, Gadfly

# const IIF = IncrementalInference
# const TU = TransformUtils

# true landmark locations
global gtl = Dict{Symbol, Vector{Float64}}()
gtl[:l1] = [1.0;-1.0;0]
gtl[:l2] = [2.0;-1.0;0]
gtl[:l3] = [3.0;-1.0;0]
gtl[:l4] = [4.0;-1.0;0]

gtl[:l11] = [1.0;1.0;0]
gtl[:l12] = [2.0;1.0;0]
gtl[:l13] = [3.0;1.0;0]
gtl[:l14] = [4.0;1.0;0]

# true pose locations
global gtp = Dict{Symbol,Vector{Float64}}()
gtp[:x0] = [0.0;0;0]
gtp[:x1] = [0.5;0;0]
gtp[:x2] = [1.0;0;0]
gtp[:x3] = [2.0;0;0]
gtp[:x4] = [2.5;0;0]
gtp[:x5] = [3.0;0;0]
gtp[:x6] = [3.0;0;0]


# calculate true measurements
global gtpt = Dict{Symbol, Vector{Float64}}()
for (ps,gp) in gtp, (ls,gl) in gtl
  gtpt[Symbol("$(ps)$(ls)")] = se2vee(SE2(gp) \ SE2(gl))
end

# Pose and Velocity noise model
global xposesig, yposesig = 0.5, 0.1
global thetasig = 0.3
global xvelsig, yvelsig = 0.2, 0.05



# Figure export folder
global currdirtime = now()
global imgdir = joinpath(ENV["HOME"], "Pictures", "testimgs", "$(currdirtime)")
mkdir(imgdir)



# construct the factor graph

global N=100
global fg = initfg()

addNode!(fg, :x0, DynPose2(ut=0))
addFactor!(fg, [:x0], DynPose2VelocityPrior(MvNormal(zeros(3),Matrix(Diagonal([0.01;0.01;0.001].^2))),
                                            MvNormal(zeros(2),Matrix(Diagonal([0.05;0.05].^2)))) )

IIF.doautoinit!(fg, :x0)

addNode!(fg, :l1, Pose2)
addFactor!(fg, [:x0;:l1], DynPose2Pose2(MvNormal(gtpt[:x0l1],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )

addNode!(fg, :l11, Pose2)
addFactor!(fg, [:x0;:l11], DynPose2Pose2(MvNormal(gtpt[:x0l11],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )


global pts = IIF.approxConv(fg, :x0l1f1, :l1, N=N)
# setVal!(fg, :l1, pts)

global pts = IIF.approxConv(fg, :x0l1f1, :x0, N=N)


global tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree, N=N)

global psid = 0
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)


# pose :x1
addNode!(fg, :x1, DynPose2(ut=1000_000))
addFactor!(fg, [:x0;:x1],VelPose2VelPose2(MvNormal(zeros(3),Matrix(Diagonal([xposesig;yposesig;thetasig].^2))),
                                          MvNormal(zeros(2),Matrix(Diagonal([xvelsig;yvelsig].^2))) ))


#
# getVal(fg, :x0)
# getVal(fg, :x1)
# getSample(IIF.getfnctype(getVert(fg, :x0x1f1, nt=:fnc)),2)
# IIF.approxConv(fg, :x0x1f1, :x1)


addFactor!(fg, [:x1;:l1], DynPose2Pose2(MvNormal(gtpt[:x1l1],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addFactor!(fg, [:x1;:l11], DynPose2Pose2(MvNormal(gtpt[:x1l11],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addNode!(fg, :l2, Pose2)
addFactor!(fg, [:x1;:l2], DynPose2Pose2(MvNormal(gtpt[:x1l2],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addNode!(fg, :l12, Pose2)
addFactor!(fg, [:x1;:l12], DynPose2Pose2(MvNormal(gtpt[:x1l12],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))

# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree, N=N)


global psid = 1
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)





# pose :x2
addNode!(fg, :x2, DynPose2(ut=2000_000))
addFactor!(fg, [:x1;:x2], VelPose2VelPose2(MvNormal(zeros(3),Matrix(Diagonal([xposesig;yposesig;thetasig].^2))),
                                           MvNormal(zeros(2),Matrix(Diagonal([xvelsig;yvelsig].^2))) ))

addFactor!(fg, [:x2;:l2], DynPose2Pose2(MvNormal(gtpt[:x2l2],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addFactor!(fg, [:x2;:l12], DynPose2Pose2(MvNormal(gtpt[:x2l12],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addNode!(fg, :l3, Pose2)
addFactor!(fg, [:x2;:l3], DynPose2Pose2(MvNormal(gtpt[:x2l3],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))

# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 2
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)






# gtp[:x3] = [2.0;0;0]
addNode!(fg, :x3, DynPose2(ut=3000_000))
addFactor!(fg, [:x2;:x3], VelPose2VelPose2(MvNormal(zeros(3),Matrix(Diagonal([xposesig;yposesig;thetasig].^2))),
                                           MvNormal(zeros(2),Matrix(Diagonal([xvelsig;yvelsig].^2))) ))

addFactor!(fg, [:x3;:l3], DynPose2Pose2(MvNormal(gtpt[:x3l3],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addNode!(fg, :l13, Pose2)
addFactor!(fg, [:x3;:l13], DynPose2Pose2(MvNormal(gtpt[:x3l13],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))

# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 3
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)






# pose :x4
# gtp[:x4] = [2.5;0;0]
addNode!(fg, :x4, DynPose2(ut=4000_000))
addFactor!(fg, [:x3;:x4], VelPose2VelPose2(MvNormal(zeros(3),Matrix(Diagonal([xposesig;yposesig;thetasig].^2))),
                                           MvNormal(zeros(2),Matrix(Diagonal([xvelsig;yvelsig].^2))) ))


addFactor!(fg, [:x4;:l3], DynPose2Pose2(MvNormal(gtpt[:x4l3],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addFactor!(fg, [:x4;:l13], DynPose2Pose2(MvNormal(gtpt[:x4l13],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addNode!(fg, :l4, Pose2)
addFactor!(fg, [:x4;:l4], DynPose2Pose2(MvNormal(gtpt[:x4l4],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addNode!(fg, :l14, Pose2)
addFactor!(fg, [:x4;:l14], DynPose2Pose2(MvNormal(gtpt[:x4l14],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))



# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 4
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)








# pose :x5
# gtp[:x5] = [3.0;0;0]
addNode!(fg, :x5, DynPose2(ut=5000_000))
addFactor!(fg, [:x4;:x5], VelPose2VelPose2(MvNormal(zeros(3),Matrix(Diagonal([xposesig;yposesig;thetasig].^2))),
                                           MvNormal(zeros(2),Matrix(Diagonal([xvelsig;yvelsig].^2))) ))

addFactor!(fg, [:x5;:l4], DynPose2Pose2(MvNormal(gtpt[:x5l4],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
addFactor!(fg, [:x5;:l14], DynPose2Pose2(MvNormal(gtpt[:x5l14],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))


# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 5
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)







# pose :x6
# gtp[:x6] = [3.0;0;0]
addNode!(fg, :x6, DynPose2(ut=6000_000))
addFactor!(fg, [:x5;:x6], VelPose2VelPose2(MvNormal(zeros(3),Matrix(Diagonal([xposesig;yposesig;thetasig].^2))),
                                           MvNormal(zeros(2),Matrix(Diagonal([xvelsig;yvelsig].^2))) ))

addFactor!(fg, [:x6;:l4], DynPose2Pose2(MvNormal(gtpt[:x6l4],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))
# addFactor!(fg, [:x6;:l14], DynPose2Pose2(MvNormal(gtpt[:x6l14],Matrix(Diagonal([0.1;0.1;0.01].^2))) ))


# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 6
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)










global coord = Coord.Cartesian(xmin=-2.0,xmax = 2.0)

global pl = plotPose2Vels(fg, :x0, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x0.png"),20cm, 10cm),pl)
global pl = plotPose2Vels(fg, :x1, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x1.png"),20cm, 10cm),pl)
global pl = plotPose2Vels(fg, :x2, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x2.png"),20cm, 10cm),pl)
global pl = plotPose2Vels(fg, :x3, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x3.png"),20cm, 10cm),pl)
global pl = plotPose2Vels(fg, :x4, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x4.png"),20cm, 10cm),pl)
global pl = plotPose2Vels(fg, :x5, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x5.png"),20cm, 10cm),pl)
global pl = plotPose2Vels(fg, :x6, coord=coord)
Gadfly.draw(PNG(joinpath(imgdir,"vel_x6.png"),20cm, 10cm),pl)


#
