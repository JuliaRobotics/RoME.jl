# additional generator tests

using Test
using RoME
using LieGroups

##

@testset "test generateGraph_TwoPoseOdo" begin
    ##

    fg = generateGraph_TwoPoseOdo()

    solveGraph!(fg)

    M = getManifold(getStateKind(fg, :x0))

    @show val1 = IIF.calcMeanMaxSuggested(fg, :x0, :simulated).suggested
    @test isapprox(M, LieGroups.identity_element(M), val1, atol = 1)
    @show val2 = IIF.calcMeanMaxSuggested(fg, :x0).suggested
    @test isapprox(M, LieGroups.identity_element(M), val2, atol = 1)

    ##
end

#
