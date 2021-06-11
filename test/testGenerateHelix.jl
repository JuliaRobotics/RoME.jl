# 

using Test
using RoME



##


@testset "test canonical helix generators" begin

##

tmp = RoME._calcHelix_T(0, 3, 25, radius=5, x_t=t->(1/3)*t)

##

fg = generateCanonicalFG_Helix2DSlew!(46, slew_x=3/2, posesperturn=15, radius=10, useMsgLikelihoods=false, Qd=diagm( [0.1;0.1;0.05].^2 ))

## # test slew in x

lastpose = sortDFG(ls(fg))[end]
@test isapprox( getPPE(fg, lastpose, :simulated).suggested , [20,0,1.465088], atol=0.001 )


##

fg = RoME.generateCanonicalFG_Helix2DSpiral!(200, graphinit=false, rate_r=0.6, rate_a=6)


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
