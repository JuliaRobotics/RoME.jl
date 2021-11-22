# 

using Test
using RoME


##

@testset "test Boxes2D generator" begin

##


fg = generateGraph_Boxes2D!(8, postpose_cb=(g,l)->@show l)

##

@test isapprox( [0;0], getPPE(fg, :x0, :simulated).suggested, atol=1e-8)
@test isapprox( [15;0], getPPE(fg, :x1, :simulated).suggested, atol=1e-8)
@test isapprox( [15;15], getPPE(fg, :x2, :simulated).suggested, atol=1e-8)
@test isapprox( [5;15], getPPE(fg, :x3, :simulated).suggested, atol=1e-8)
@test isapprox( [5;0], getPPE(fg, :x4, :simulated).suggested, atol=1e-8)
@test isapprox( [20;0], getPPE(fg, :x5, :simulated).suggested, atol=1e-8)
@test isapprox( [20;15], getPPE(fg, :x6, :simulated).suggested, atol=1e-8)
@test isapprox( [10;15], getPPE(fg, :x7, :simulated).suggested, atol=1e-8)
@test isapprox( [10;0], getPPE(fg, :x8, :simulated).suggested, atol=1e-8)

##

solveTree!(fg)

##

@test isapprox( [0;0], getPPE(fg, :x0, :simulated).suggested, atol=1e-8)
@test isapprox( [15;0], getPPE(fg, :x1, :simulated).suggested, atol=1e-8)
@test isapprox( [15;15], getPPE(fg, :x2, :simulated).suggested, atol=1e-8)
@test isapprox( [5;15], getPPE(fg, :x3, :simulated).suggested, atol=1e-8)
@test isapprox( [5;0], getPPE(fg, :x4, :simulated).suggested, atol=1e-8)
@test isapprox( [20;0], getPPE(fg, :x5, :simulated).suggested, atol=1e-8)
@test isapprox( [20;15], getPPE(fg, :x6, :simulated).suggested, atol=1e-8)
@test isapprox( [10;15], getPPE(fg, :x7, :simulated).suggested, atol=1e-8)
@test isapprox( [10;0], getPPE(fg, :x8, :simulated).suggested, atol=1e-8)

##

end

@testset "test canonical helix generators" begin

##

# moved to IIF
# tmp = calcHelix_T(0, 3, 25, radius=5, xr_t=t->(1/3)*t)

##

cb(fg_, lp) = @show lp, length(ls(fg_))

fg = generateGraph_Helix2DSlew!(46, slew_x=2/3, posesperturn=15, radius=10, Qd=diagm( [0.1;0.1;0.05].^2 ), postpose_cb=cb, solverParams=SolverParams(useMsgLikelihoods=false))

## # test slew in x

lastpose = sortDFG(ls(fg))[end]
@test isapprox( getPPE(fg, lastpose, :simulated).suggested , [20,0,1.465088], atol=0.001 )


##

fg = RoME.generateGraph_Helix2DSpiral!(200, rate_r=0.6, rate_a=6, radius=100, solverParams=SolverParams(graphinit=false))


##
end


@testset "test extending fg with helix generator" begin
##

fg = generateGraph_Helix2D!(5, posesperturn=15, radius=10, solverParams=SolverParams(graphinit=false))

@test !getSolverParams(fg).graphinit

# check the helix is being constructed in a consistent way, using ppe solveKey :simulated
vars = sortDFG(ls(fg))
ppes = [
  [0.0, 0.0, 1.5707963267948966],
  [0.8645454235739924, 4.067366430758004, 1.151917276019672],
  [3.3086939364114176, 7.431448254773942, 0.7330382545911657],
  [6.909830056250526, 9.510565162951536, 0.31415923447063226],
  [11.045284632676536, 9.945218953682733, -0.10471978645923721],
]

for (i,v) in enumerate(vars)
  @test isapprox(ppes[i], getPPE(fg, v, :simulated).suggested; atol=1e-5 ) 
end


# check that the graph can be expanded with the same generator function

generateGraph_Helix2D!(5, dfg=fg, posesperturn=15, radius=10, solverParams=SolverParams(graphinit=false))
# check that no new variables were added
@test length(ls(fg)) == length(vars)
@test !exists(fg, :x5)

##

generateGraph_Helix2D!(6, dfg=fg, posesperturn=15, radius=10, solverParams=SolverParams(graphinit=false))
# check that only one variable has been added
@test length(ls(fg)) == length(vars) + 1 
@test exists(fg, :x5)

# check the simulated location of the new variable is also at the right location
push!(vars, :x5)
push!(ppes, [ 15.0;8.660254037844387;-0.5235988055902416] )

@test isapprox(ppes[end], getPPE(fg, :x5, :simulated).suggested; atol=1e-5 )


##
end


##

# using TensorCast
# using Gadfly
# Gadfly.set_default_plot_size(35cm,25cm)

# ##

# vals = getPPE.(fg, (ls(fg) |> sortDFG), :simulated) .|> x-> x.suggested
# @cast XYT[d,p] := vals[p][d]

# Gadfly.plot(x=XYT[1,:],y=XYT[2,:], Geom.path)
  
##
