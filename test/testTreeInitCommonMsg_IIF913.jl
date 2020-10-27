
using Test
using Statistics
using RoME

@testset "Testing tree init prior usage" begin

fg = generateCanonicalFG_Circle(4;graphinit=true)
# deleteVariable!.(fg, [:x3, :x4, :l1])


# fg.solverParams.showtree = true
# fg.solverParams.drawtree = true
# fg.solverParams.dbg = true
fg.solverParams.useMsgLikelihoods = true

# now move the graph prior somewhere else (ie all init is wrong)
deleteFactor!(fg, :x0f1)
prpo = PriorPose2(MvNormal([5.,0.,0.0], 0.01*Matrix{Float64}(LinearAlgebra.I,3,3)))
addFactor!(fg, [:x0], prpo)

solveTree!(fg);

X4 = getBelief(fg, :x4) |> getPoints


@error("Must first fix IIF #913")
# @test_broken 2.0 < Statistics.mean(X4[1,:])
@test -1.5 < Statistics.mean(X4[2,:]) < 1.5


end