using RoME
using Statistics
using JSON2

using Test

##

@testset "basic Point2Point2 test" begin

fg = initfg()

addVariable!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), Matrix{Float64}(LinearAlgebra.I, 2,2))))

addVariable!(fg, :x1, Point2)
addFactor!(fg, [:x0;:x1], Point2Point2(MvNormal([10;0.0], Matrix{Float64}(LinearAlgebra.I, 2,2))))

tree = solveTree!(fg)
# tree = wipeBuildNewTree!(fg)
# inferOverTree!(fg, tree)

@test isapprox(mean(getVal(fg, :x0)), [0,0], atol=1.0)
@test isapprox(mean(getVal(fg, :x1)), [10,0], atol=1.0)

end

##

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


@testset "test Point2Point2Range..." begin

##

N=100 # return to 200
fg = initfg()
# getSolverParams(fg).inflation = 50.0
getSolverParams(fg).graphinit = false

addVariable!(fg, :x0, Point2, N=N)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], diagm(ones(2)) )), graphinit=false)

addVariable!(fg, :x1, Point2, N=N)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], diagm(ones(2)) )), graphinit=false)

addVariable!(fg, :l1, Point2, N=N)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)) , graphinit=false)
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)) , graphinit=false)


##

@warn("Point2Point2 range 2 mode, allow 3 attempts until IIF #1010 is completed")
TP = false
for ic in 1:3
  tree = solveTree!(fg)

  # mode 1
  @cast l1_val[j,i] := getVal(fg, :l1)[i][j]
  @show T1 = (0.05*N < sum( 90 .< l1_val[1,:] .< 110 ) && 0.05*N < sum( 90 .< l1_val[2,:] .< 110 ))
  # mode 2
  @show T2 = (0.05*N < sum( -10 .< l1_val[1,:] .< 10 ) && 0.05*N < sum( -10 .< l1_val[2,:] .< 10 ))
  TP |= T1 && T2
  TP && break
end

@test TP

# @test 0.05*N < sum( 90 .< getVal(fg, :l1)[1,:] .< 110 )
# @test 0.05*N < sum( -10 .< getVal(fg, :l1)[2,:] .< 10 )

# @test 0.05*N < sum( -10 .< getVal(fg, :l1)[1,:] .< 10 )
# @test 0.05*N < sum( 90 .< getVal(fg, :l1)[2,:] .< 110 )
@cast l1_val[j,i] := getVal(fg, :l1)[i][j]

voidsel1 =  10.0 .< l1_val[1,:]
@test sum( l1_val[2,voidsel1] .< 70 ) < 0.35*N

voidsel2 =  10.0 .< l1_val[2,:]
@test sum( l1_val[1,voidsel2] .< 70 ) < 0.35*N

@test sum( 120 .< abs.(l1_val[1,:]) ) < 0.35*N
@test sum( 120 .< abs.(l1_val[2,:]) ) < 0.35*N

##

end


##
@testset "Unpack of Point2Point2Range, #563" begin
##

fg = initfg()

addVariable!(fg, :x0, Point2)
addVariable!(fg, :l3, Point2)

point2point2rangeString = "{\"label\":\"x0l3f1\",\"_version\":\"0.18.1\",\"_variableOrderSymbols\":[\"x0\",\"l3\"],\"data\":{\"eliminated\":false,\"potentialused\":false,\"edgeIDs\":[],\"fnc\":{\"Z\":{\"_type\":\"IncrementalInference.PackedNormal\",\"mu\":89.44271909999159,\"sigma\":3.0}},\"multihypo\":[],\"certainhypo\":[1,2],\"nullhypo\":0.0,\"solveInProgress\":0,\"inflation\":5.0},\"tags\":[\"FACTOR\"],\"timestamp\":\"2022-03-26T01:24:44.373-05:00\",\"nstime\":\"0\",\"fnctype\":\"Point2Point2Range\",\"solvable\":1}"
jback = JSON2.read(point2point2rangeString, Dict{String, Any})
f = unpackFactor(fg, jback)


##
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
