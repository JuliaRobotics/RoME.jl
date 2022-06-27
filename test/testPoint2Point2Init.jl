
using RoME
using Test
using TensorCast

##

@testset "test tree based autoinit on small bimodal point2 trilateration..." begin

##

N=150
fg = initfg()
getSolverParams(fg).graphinit=false
getSolverParams(fg).treeinit=false
getSolverParams(fg).useMsgLikelihoods = true
getSolverParams(fg).gibbsIters = 5
# getSolverParams(fg).inflation = 10.0


addVariable!(fg, :x0, Point2, N=N)
addFactor!(fg, [:x0], PriorPoint2(MvNormal([100.0;0], [1.0 0; 0 1])), graphinit=false)

addVariable!(fg, :x1, Point2, N=N)
addFactor!(fg, [:x1], PriorPoint2(MvNormal([0.0;100.0], [1.0 0; 0 1])), graphinit=false)

addVariable!(fg, :l1, Point2, N=N)
addFactor!(fg, [:x0;:l1], Point2Point2Range(Normal(100.0, 1.0)), graphinit=false)
addFactor!(fg, [:x1;:l1], Point2Point2Range(Normal(100.0, 1.0)), graphinit=false)

## what would 

eo = [:x1; :x0; :l1]

## Mock clique solve

# tree = buildTreeReset!(fg, eo)
# # from clique 2
# sfg2 = buildCliqSubgraph(fg, tree, :x1)

# prior first

# pts, = predictbelief(sfg2,:x1, :)
# initVariable!(sfg2, :x1, pts)
# pts, = predictbelief(sfg2, :l1, :)
# initVariable!(sfg2, :l1, pts)


##

# plotKDE(getBelief(sfg2,:l1))

##

i = 3
while 0 < i
##
  tree = solveTree!(fg, eliminationOrder=eo);
##
  try
    # 
    @cast l1_vals[j,i] := getVal(fg, :l1)[i][j]
    
    @test 5 < sum( 80 .< l1_vals[1,:] .< 120 )
    @test 5 < sum( -20 .< l1_vals[2,:] .< 20 )
    
    @test 5 < sum( -20 .< l1_vals[1,:] .< 20 )
    @test 5 < sum( 80 .< l1_vals[2,:] .< 120 )
    
    @cast x0_vals[j,i] := getVal(fg, :x0)[i][j]
    
    @test 80 < sum( 80 .< x0_vals[1,:] .< 120 )
    @test 80 < sum( -20 .< x0_vals[2,:] .< 20 )
    
    @cast x1_vals[j,i] := getVal(fg, :x1)[i][j]
    
    @test 80 < sum( -20 .< x1_vals[1,:] .< 20 )
    @test 80 < sum( 80 .< x1_vals[2,:] .< 120 )

    break
    #
  catch
    i -= 1
    i < 1 ? error("3 attempts failed") : nothing
  end
end



## testset
end

##

# using RoMEPlotting
# Gadfly.set_default_plot_size(25cm,20cm)

##

# plotKDE(fg, [:x0;:x1]) # |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)
# plotKDE(fg, :l1) # |> SVG("/tmp/test.svg") || @async run(`eog /tmp/test.svg`)

##

# IIF._getCCW(fg, :x1l1f1).inflation = 100.0
# pts = approxConv(fg, :x1l1f1, :l1)
# initVariable!(fg, :l1, pts)

# plotKDE(manikde!(Point2, pts))

##