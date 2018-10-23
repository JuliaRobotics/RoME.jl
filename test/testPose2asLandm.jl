# build test case for using Pose2Pose2 to landmarks -- to help AprilTag

using RoME, Distributions
using RoMEPlotting, Gadfly

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


# calculate true measurements
global gtpt = Dict{Symbol, Vector{Float64}}()
for (ps,gp) in gtp, (ls,gl) in gtl
  gtpt[Symbol("$(ps)$(ls)")] = se2vee(SE2(gp) \ SE2(gl))
end




# Figure export folder
global currdirtime = now()
global imgdir = joinpath(ENV["HOME"], "Pictures", "testimgs", "$(currdirtime)")
mkdir(imgdir)



# construct the factor graph

global N=75
global fg = initfg()

addNode!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3),Matrix(Diagonal([0.01;0.01;0.001].^2)))))

addNode!(fg, :l1, Pose2)
addFactor!(fg, [:x0;:l1], Pose2Pose2(MvNormal(gtpt[:x0l1],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l11, Pose2)
addFactor!(fg, [:x0;:l11], Pose2Pose2(MvNormal(gtpt[:x0l11],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)

global psid = 0
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)



# pose :x1
addNode!(fg, :x1, Pose2)
addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal(zeros(3),Matrix(Diagonal([0.4;0.1;0.4].^2)))) )

addFactor!(fg, [:x1;:l1], Pose2Pose2(MvNormal(gtpt[:x1l1],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addFactor!(fg, [:x1;:l11], Pose2Pose2(MvNormal(gtpt[:x1l11],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l2, Pose2)
addFactor!(fg, [:x1;:l2], Pose2Pose2(MvNormal(gtpt[:x1l2],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l12, Pose2)
addFactor!(fg, [:x1;:l12], Pose2Pose2(MvNormal(gtpt[:x1l12],Matrix(Diagonal([0.1;0.1;0.01].^2)))) )

# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 1
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)





# pose :x2
addNode!(fg, :x2, Pose2)
addFactor!(fg, [:x1;:x2], Pose2Pose2(MvNormal(zeros(3), Matrix(Diagonal([0.4;0.1;0.4].^2)))) )

addFactor!(fg, [:x2;:l2], Pose2Pose2(MvNormal(gtpt[:x2l2], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addFactor!(fg, [:x2;:l12], Pose2Pose2(MvNormal(gtpt[:x2l12], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l3, Pose2)
addFactor!(fg, [:x2;:l3], Pose2Pose2(MvNormal(gtpt[:x2l3], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )

# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 2
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)






# gtp[:x3] = [2.0;0;0]
addNode!(fg, :x3, Pose2)
addFactor!(fg, [:x2;:x3], Pose2Pose2(MvNormal(zeros(3), Matrix(Diagonal([0.4;0.1;0.4].^2)))) )

addFactor!(fg, [:x3;:l3], Pose2Pose2(MvNormal(gtpt[:x3l3], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l13, Pose2)
addFactor!(fg, [:x3;:l13], Pose2Pose2(MvNormal(gtpt[:x3l13], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )

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
addNode!(fg, :x4, Pose2)
addFactor!(fg, [:x3;:x4], Pose2Pose2(MvNormal(zeros(3), Matrix(Diagonal([0.4;0.1;0.4].^2)))) )


addFactor!(fg, [:x4;:l3], Pose2Pose2(MvNormal(gtpt[:x4l3], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addFactor!(fg, [:x4;:l13], Pose2Pose2(MvNormal(gtpt[:x4l13], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l4, Pose2)
addFactor!(fg, [:x4;:l4], Pose2Pose2(MvNormal(gtpt[:x4l4], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addNode!(fg, :l14, Pose2)
addFactor!(fg, [:x4;:l14], Pose2Pose2(MvNormal(gtpt[:x4l14], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )



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
addNode!(fg, :x5, Pose2)
addFactor!(fg, [:x4;:x5], Pose2Pose2(MvNormal(zeros(3), Matrix(Diagonal([0.4;0.1;0.4].^2)))) )

addFactor!(fg, [:x5;:l4], Pose2Pose2(MvNormal(gtpt[:x5l4], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )
addFactor!(fg, [:x5;:l14], Pose2Pose2(MvNormal(gtpt[:x5l14], Matrix(Diagonal([0.1;0.1;0.01].^2)))) )


# writeGraphPdf(fg)

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)


global psid = 5
global pl = drawPosesLandms(fg, spscale=0.1, drawhist=false)#,   meanmax=:mean,xmin=-3,xmax=6,ymin=-5,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"x$(psid).png"),30cm, 25cm),pl)
global pl = drawPosesLandms(fg, spscale=0.1)#,   meanmax=:mean,xmin=-3,xmax=3,ymin=-2,ymax=2);
Gadfly.draw(PNG(joinpath(imgdir,"hist_x$(psid).png"),30cm, 25cm),pl)




#
