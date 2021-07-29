# Basic tests of poses at zero

using RoME
using Test
using Manifolds
using TensorCast

@testset "basic pose2 trivial case without forcing autoinit..." begin

fg = initfg()

cov = Matrix(Diagonal(1e-4.*[1;1;1]))

# first pose position
addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3),cov)))

# second Pose position with odo
addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal(zeros(3),cov)))

# third pose position with odo
addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x1;:x2], Pose2Pose2(MvNormal(zeros(3),cov)))

#

Xc_badval = 0.01.*randn(3,100)
Xc_badval[1,:] .-= 5.0
Xc_badval[2,:] .-= 2.0
Xc_badval[3,:] .+= 0.5

badval = map((Xc)->DFG.getPoint(Pose2, Xc), eachcol(Xc_badval))

setValKDE!(fg, :x2, manikde!(badval, Pose2))

N = 100
# batchSolve!(fg, N=N)
tree, smt, hist = solveTree!(fg)

for xx in [:x0; :x1; :x2]
  pts = getVal(fg, xx)
  pts_μ = mean(getManifold(Pose2), pts)
  isapprox(getManifold(Pose2), pts_μ, ProductRepr([0.,0], [1. 0; 0 1]), atol=0.01)


  Xpts = getCoordinates.(Ref(Pose2), pts) 
  @cast Xptsarr[j,i] := Xpts[i][j]

  @test 0.95*N < sum( -1.0 .< Xptsarr[1,:] .< 1.0 )
  @test 0.95*N < sum( -1.0 .< Xptsarr[2,:] .< 1.0 )
  @test 0.95*N < sum( -0.5 .< Xptsarr[3,:] .< 0.5 )
end


end







@testset "basic pose2 with forcing bad initialization..." begin

fg = initfg()

cov = Matrix(Diagonal(1e-4.*[1;1;1]))

# first pose position
addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3),cov)))

# second Pose position with odo
addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal(zeros(3),cov)))

# third pose position with odo
addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x1;:x2], Pose2Pose2(MvNormal(zeros(3),cov)))

ensureAllInitialized!(fg)


Xc_badval = 0.000001.*randn(3,100)
Xc_badval[1,:] .-= 5.0
Xc_badval[2,:] .-= 2.0
Xc_badval[3,:] .+= 0.5
badval = map((Xc)->DFG.getPoint(Pose2, Xc), eachcol(Xc_badval))

setValKDE!(fg, :x2, manikde!(badval, Pose2))

# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)

N = 100
getSolverParams(fg).N = N
solveTree!(fg)



for xx in [:x0; :x1; :x2]
  manipts = getVal(fg, xx)
  manicrd = getCoordinates.(Ref(Pose2), manipts) 
  @cast pts[j,i] := manicrd[i][j]
  
  @test 0.95*N < sum( -0.5 .< pts[1,:] .< 0.5 )
  @test 0.95*N < sum( -1.0 .< pts[2,:] .< 1.0 )
  @test 0.95*N < sum( -1.0 .< pts[3,:] .< 1.0 )
end

end



# using RoMEPlotting
#
# drawPoses(fg, spscale=0.5)
#
# spyCliqMat(tree, :x2)



@testset "test basic banana (split)..." begin


fg = initfg()

addVariable!(fg, :x0, Pose2())
addVariable!(fg, :x1, Pose2())

addFactor!(fg,
           [:x0],
           PriorPose2(MvNormal(zeros(3), Matrix(Diagonal([0.01;0.01;1.0].^2)))),
           graphinit=false )
#

addFactor!(fg,
           [:x0;:x1],
           Pose2Pose2(MvNormal([1.0;0;0], Matrix(Diagonal([0.01;0.01;0.01].^2)))),
           graphinit=false )
#

solveTree!(fg)


manipts = getPoints(getKDE(fg, :x0))
manicrd = getCoordinates.(Ref(Pose2), manipts) 
@cast pts[j,i] := manicrd[i][j]

N = size(pts,2)

@test 0.7*N < sum(abs.(pts[1,:]) .< 0.1)
@test 0.7*N < sum(abs.(pts[2,:]) .< 0.1)
@test 0.7*N < sum(abs.(pts[3,:]) .< 2.0)



manipts = getPoints(getKDE(fg, :x1))
manicrd = getCoordinates.(Ref(Pose2), manipts) 
@cast pts[j,i] := manicrd[i][j]

@test 0.7*N < sum(0.0 .< pts[1,:])
@test 0.7*N < sum( 0.5 .< sqrt.(sum(pts[1:2,:].^2, dims=1)) .< 1.5)
@test 0.7*N < sum(abs.(pts[3,:]) .< 2.0)



#
# using RoMEPlotting
# using Gadfly
# Gadfly.set_default_plot_size(35cm,25cm)
#
# plotKDE(fg, ls(fg))


end
