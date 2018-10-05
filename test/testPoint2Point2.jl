using RoME
using Statistics
# , Distributions

using Test

@testset "basic Point2Point2 test" begin

global fg = initfg()

addNode!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), Matrix{Float64}(LinearAlgebra.I, 2,2))))

addNode!(fg, :x1, Point2)
addFactor!(fg, [:x0;:x1], Point2Point2(MvNormal([10;0.0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)

@test sum( abs.(Statistics.mean(getVal(fg, :x0),dims=2) .- [0.0;0]) .< [0.5;0.5]) == 2
@test sum( abs.(Statistics.mean(getVal(fg, :x1),dims=2) .- [10.0;0]) .< [0.5;0.5]) == 2

end


@testset "test Point2Point2Range{T}..." begin

global N=200
global fg = initfg()

addNode!(fg, :x0, Point2, N=N)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

addNode!(fg, :x1, Point2, N=N)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

addNode!(fg, :l1, Point2, N=N)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)))
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)))


global tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree, N=N)

@test sum( 90 .< getVal(fg, :l1)[1,:] .< 110 ) > 32
@test sum( -10 .< getVal(fg, :l1)[2,:] .< 10 ) > 32

@test sum( -10 .< getVal(fg, :l1)[1,:] .< 10 ) > 32
@test sum( 90 .< getVal(fg, :l1)[2,:] .< 110 ) > 32

global voidsel1 =  10.0 .< getVal(fg, :l1)[1,:]
@test sum( getVal(fg, :l1)[2,voidsel1] .< 80 ) < 15

global voidsel2 =  10.0 .< getVal(fg, :l1)[2,:]
@test sum( getVal(fg, :l1)[1,voidsel2] .< 80 ) < 15

@test sum( 120 .< abs.(getVal(fg, :l1)[1,:]) ) < 15
@test sum( 120 .< abs.(getVal(fg, :l1)[2,:]) ) < 15


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
