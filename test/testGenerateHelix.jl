# 

using Test
using RoME


##

@testset "test Boxes2D generator" begin

##


fg = generateCanonicalFG_Boxes2D!(8, postpose_cb=(g,l)->@show l)

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

tmp = calcHelix_T(0, 3, 25, radius=5, xr_t=t->(1/3)*t)

##

cb(fg_, lp) = @show lp, length(ls(fg_))

fg = generateCanonicalFG_Helix2DSlew!(46, slew_x=2/3, posesperturn=15, radius=10, useMsgLikelihoods=false, Qd=diagm( [0.1;0.1;0.05].^2 ), postpose_cb=cb)

## # test slew in x

lastpose = sortDFG(ls(fg))[end]
@test isapprox( getPPE(fg, lastpose, :simulated).suggested , [20,0,1.465088], atol=0.001 )


##

fg = RoME.generateCanonicalFG_Helix2DSpiral!(200, graphinit=false, rate_r=0.6, rate_a=6, radius=100)


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
