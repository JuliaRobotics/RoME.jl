
using RoME
using RoMEPlotting


fg = initfg()

addVariable!(fg, :l2, Point2)
addFactor!(fg, [:l2;], Prior(MvNormal([0.0;0],Matrix(Diagonal([0.01;0.01])))), autoinit=false)

addVariable!(fg, :x13, Pose2)
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x13;:l2], p2br, autoinit=false)




tree = batchSolve!(fg, treeinit=true, drawpdf=true)


drawPosesLandms(fg)

plotPose(fg, :x13)



#
