# 

using Test
using RoME



##


@testset "test canonical helix generators" begin

##

tmp = RoME._calcHelix_T(0, 3, 25, radius=5, x_t=t->(1/3)*t)

##

fg = RoME.generateCanonicalFG_Helix2D!(46, posesperturn=15, radius=10, useMsgLikelihoods=false, Qd=diagm( [0.1;0.1;0.05].^2 ))

##

lastpose = sortDFG(ls(fg))[end]

@test isapprox( getPPE(fg, lastpose, :simulated).suggested , [10,0,1.517794], atol=0.001 )

##

end

##

# using TensorCast
# using Gadfly

# ##

# vals = getPPE.(fg, (ls(fg) |> sortDFG), :simulated) .|> x-> x.suggested
# @cast XYT[d,p] := vals[p][d]

# Gadfly.plot(x=XYT[1,:],y=XYT[2,:], Geom.path)
  
##
