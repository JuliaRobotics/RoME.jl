# 

using Test
using RoME



##


@testset "test canonical helix generators" begin


st = RoME._calcHelix2DTurnsX(5, radius=10, N_ppt=20, runback=5/7)

fg = RoME.generateCanonicalFG_Helix2D!(28, useMsgLikelihoods=false, Qd=diagm( [0.1;0.1;0.05].^2 ))


end


