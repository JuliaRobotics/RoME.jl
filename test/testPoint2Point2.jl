using RoME
using Statistics
# , Distributions

using Test

@testset "basic Point2Point2 test" begin

global fg = initfg()

addVariable!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), Matrix{Float64}(LinearAlgebra.I, 2,2))))

addVariable!(fg, :x1, Point2)
addFactor!(fg, [:x0;:x1], Point2Point2(MvNormal([10;0.0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

tree, smt, hist = solveTree!(fg)
# global tree = wipeBuildNewTree!(fg)
# inferOverTree!(fg, tree)

@test sum( abs.(Statistics.mean(getVal(fg, :x0),dims=2) .- [0.0;0]) .< [0.5;0.5]) == 2
@test sum( abs.(Statistics.mean(getVal(fg, :x1),dims=2) .- [10.0;0]) .< [0.5;0.5]) == 2

end


# drawGraph(fg, show=true)
#
# using RoMEPlotting, Gadfly, Fontconfig, Cairo
#
# drawTree(tree, show=true)
#
# stuff = treeProductUp(fg, tree, :x0, :l1)
# stuff = treeProductUp(fg, tree, :l1, :l1)
#
# plotKDE(kde!(stuff[1]), levels=3)
# # drawLandms(fg) #, regexLandmark=r"x") |> PDF("/tmp/test.pdf")
# plotKDE(fg, ls(fg))


@testset "test Point2Point2Range{T}..." begin

global N=100 # return to 200
global fg = initfg()

addVariable!(fg, :x0, Point2, N=N)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

addVariable!(fg, :x1, Point2, N=N)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

addVariable!(fg, :l1, Point2, N=N)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)))
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)))


tree, smt, hist = solveTree!(fg)
# global tree = wipeBuildNewTree!(fg)
# inferOverTree!(fg, tree, N=N)


@test 0.15*N < sum( 90 .< getVal(fg, :l1)[1,:] .< 110 )
@test 0.15*N < sum( -10 .< getVal(fg, :l1)[2,:] .< 10 )

@test 0.15*N < sum( -10 .< getVal(fg, :l1)[1,:] .< 10 )
@test 0.15*N < sum( 90 .< getVal(fg, :l1)[2,:] .< 110 )

global voidsel1 =  10.0 .< getVal(fg, :l1)[1,:]
@test sum( getVal(fg, :l1)[2,voidsel1] .< 70 ) < 0.3*N

global voidsel2 =  10.0 .< getVal(fg, :l1)[2,:]
@test sum( getVal(fg, :l1)[1,voidsel2] .< 70 ) < 0.3*N

@test sum( 120 .< abs.(getVal(fg, :l1)[1,:]) ) < 0.3*N
@test sum( 120 .< abs.(getVal(fg, :l1)[2,:]) ) < 0.3*N


end

#
#
# # debugging
#

# setVal!(fg, :x0, zeros(2,1))
# setVal!(fg, :l1, zeros(2,1))
# #
# using RoMEPlotting, KernelDensityEstimatePlotting
#
# plotKDE(fg, :x0)
#
#
# pts1 = IIF.approxConv(fg, :x0l1f1, :l1, N=N)
#
# plotKDE(KDE.kde!(pts1))
#
#
# pts2 = IIF.approxConv(fg, :x1l1f1, :l1, N=N)
#
# plotKDE([KDE.kde!(pts1); KDE.kde!(pts2)])




# #
# #
# #
#
# #



# using RoMEPlotting, KernelDensityEstimatePlotting
#
# plotKDE(fg, [:x0; :x1; :l1], levels=3)
#
# stuff = IIF.localProduct(fg, :l1)
# stuff
#
# plotKDE(stuff[2], levels=3, c=["red";"cyan"])








##
