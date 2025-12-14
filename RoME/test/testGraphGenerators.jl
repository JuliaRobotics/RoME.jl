# additional generator tests

using Test
using RoME

##

@testset "test generateGraph_TwoPoseOdo" begin
    ##

    fg = generateGraph_TwoPoseOdo()

    solveGraph!(fg)

    @test isapprox([0; 0; 0.0], IIF.calcMeanMaxSuggested(fg, :x0, :simulated).suggested, atol = 1)
    @test isapprox([0; 0; 0.0], IIF.calcMeanMaxSuggested(fg, :x0).suggested, atol = 1)

    ##
end

#
