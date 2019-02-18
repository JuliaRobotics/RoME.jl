# Basic tests of poses at zero

using RoME
using Test

@testset "basic that pose2 pushed back to zero in trivial case solve..." begin

fg = initfg()

cov = Matrix(Diagonal(1e-4.*[1;1;1]))

# first pose position
addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0], Prior(MvNormal(zeros(3),cov)))

# second Pose position with odo
addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal(zeros(3),cov)))

# third pose position with odo
addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x1;:x2], Pose2Pose2(MvNormal(zeros(3),cov)))

#

badval = 0.01.*randn(3,100)
badval[1,:] .-= 5.0
badval[2,:] .-= 2.0
badval[3,:] .+= 0.5
setValKDE!(fg, :x2, kde!(badval))

N = 100
batchSolve!(fg, N=N)

pts = getVal(fg, :x2)

@test 0.95*N < sum( -0.5 .< pts[1,:] .< 0.5 )
@test 0.95*N < sum( -0.5 .< pts[2,:] .< 0.5 )
@test 0.95*N < sum( -0.5 .< pts[3,:] .< 0.5 )


end


# using RoMEPlotting
#
# drawPoses(fg, spscale=0.5)






#
