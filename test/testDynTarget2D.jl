using Caesar

fg = initfg()

cvf = VelPoint2VelPoint2(MvNormal([0.0;0;0;0], Matrix(Diagonal([7.0;7.0;0.1;0.1].^2))))
# cvf = DynPoint2DynPoint2(MvNormal([0.0;0;0;0], Matrix(Diagonal([7.0;7.0;0.1;0.1].^2))))
br = Pose2Point2BearingRange(Normal(pi/2, 0.05), Normal(50.0, 1.0))

addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x1], PriorPose2(MvNormal([0.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l1, DynPoint2(ut=0))
addFactor!(fg, [:x1;:l1], br)

addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x2], PriorPose2(MvNormal([10.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l2, DynPoint2(ut=1_000_000))
addFactor!(fg, [:l1;:l2], cvf)
addFactor!(fg, [:x2;:l2], br)

addVariable!(fg, :x3, Pose2)
addFactor!(fg, [:x3], PriorPose2(MvNormal([20.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l3, DynPoint2(ut=2_000_000))
addFactor!(fg, [:l2;:l3], cvf)
addFactor!(fg, [:x3;:l3], br)

addVariable!(fg, :x4, Pose2)
addFactor!(fg, [:x4], PriorPose2(MvNormal([30.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l4, DynPoint2(ut=3_000_000))
addFactor!(fg, [:l3;:l4], cvf)
addFactor!(fg, [:x4;:l4], br)


# specific parameters
getSolverParams(fg).drawtree = true

# and solve
tree, smt, hist = solveTree!(fg)



## and plot
using RoMEPlotting
Gadfly.set_default_plot_size(40cm,25cm)



plotKDE(fg, sort(union(ls(fg,r"x"),ls(fg,r"l"))), dims=[1;2], title="Positions")

plotKDE(fg, sort(ls(fg,r"l")), dims=[3;4], title="velocities")










fg = initfg()

cvf = DynPoint2DynPoint2(MvNormal([0.0;0;0;0], Matrix(Diagonal([7.0;7.0;0.1;0.1].^2))))
br = Pose2Point2BearingRange(Normal(-pi/4, 0.05), Normal(20.0, 1.0))

addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x1], PriorPose2(MvNormal([0.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l1, DynPoint2(ut=0))
addFactor!(fg, [:x1;:l1], br)

addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x2], PriorPose2(MvNormal([10.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l2, DynPoint2(ut=1_000_000))
addFactor!(fg, [:l1;:l2], cvf)
addFactor!(fg, [:x2;:l2], br)

addVariable!(fg, :x3, Pose2)
addFactor!(fg, [:x3], PriorPose2(MvNormal([20.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l3, DynPoint2(ut=2_000_000))
addFactor!(fg, [:l2;:l3], cvf)
addFactor!(fg, [:x3;:l3], br)

addVariable!(fg, :x4, Pose2)
addFactor!(fg, [:x4], PriorPose2(MvNormal([30.0;0;0],Matrix(Diagonal([0.1;0.1;0.05].^2)))))
addVariable!(fg, :l4, DynPoint2(ut=3_000_000))
addFactor!(fg, [:l3;:l4], cvf)
addFactor!(fg, [:x4;:l4], br)


# specific parameters
getSolverParams(fg).drawtree = true

# and solve
tree, smt, hist = solveTree!(fg)



## and plot
using RoMEPlotting
Gadfly.set_default_plot_size(40cm,25cm)



plotKDE(fg, sort(union(ls(fg,r"x"),ls(fg,r"l"))), dims=[1;2], title="Positions")

plotKDE(fg, sort(ls(fg,r"l")), dims=[3;4], title="velocities")













#
#
# ## Test VelPoint2VelPoint2 factor
#
# res = [0.0;]
#
# userdata = FactorMetadata()
# userdata.variableuserdata = (DynPoint2(ut=0),DynPoint2(ut=1000000))
# idx = 1
# # meas = getSample(cvf)
# meas = (reshape([10;0;0;0.0],4,1),)
#
# Xi = zeros(4,1)
# Xj = zeros(4,1)
# Xj[1,1] = 10.0
#
# cost = (vx,vy) -> cvf(res,userdata,idx,meas,reshape([vy;Xi[2,1];vx;Xi[4,1]],4,1), Xj )
#
#
# cost(0.0,0.0)
# cost(5.0,0.0)
# cost(10.0,0.0)
#
#
#
# ## and look at the cost function
#
# using Gadfly
#
# vx=-20:0.1:20
# vy=-20:0.1:20
#
# Gadfly.plot(x=vx,y=vy,z=cost, Geom.contour)
#
#
#
# 0
#
# # new factor

#
#
# """
#     $(TYPEDEF)
#
# Partial prior belief on Z, Roll, and Pitch of a `Pose3`.
# """
# mutable struct PriorPose3ZRP{T1,T2} <: IncrementalInference.FunctorSingleton where {T1 <: SamplableBelief, T2 <: SamplableBelief}
#   z::T1
#   rp::T2
#   partial::Tuple
#   PriorPose3ZRP{T1,T2}() where {T1, T2} = new{T1,T2}()
#   PriorPose3ZRP{T1,T2}(z::T1,rp::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(z, rp, (3,4,5))
# end
# PriorPose3ZRP(z::T1,rp::T2) where {T1 <: SamplableBelief, T2 <: SamplableBelief} = PriorPose3ZRP{T1,T2}(z, rp)
# function getSample(pprz::PriorPose3ZRP, N::Int=1)
#   return ([rand(pprz.z,N)[:]';rand(pprz.rp,N)], )
# end
#
#



#
