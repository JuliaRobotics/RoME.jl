# beehive larger test

using Test
using RoME
using Distributed

while nprocs() < 3
  addprocs(1)
end

using RoME
@everywhere using RoME

##

@testset "sequential beehive, repeated growing and solving" begin

##

fg = RoME.generateCanonicalFG_Honeycomb!(7, graphinit=true, useMsgLikelihoods = true)
tree, _, _ = solveTree!(fg);

fg = RoME.generateCanonicalFG_Honeycomb!(14, graphinit=true, fg=fg)
tree, _, _ = solveTree!(fg, tree);

tree_ = deepcopy(tree); fg_ = deepcopy(fg);
fg = RoME.generateCanonicalFG_Honeycomb!(21, graphinit=true, fg=fg)
tree, _, _ = solveTree!(fg  , tree);

# fg = RoME.generateCanonicalFG_Honeycomb!(28, graphinit=true, fg=fg)
# tree, _, _ = solveTree!(fg, tree);

# fg = RoME.generateCanonicalFG_Honeycomb!(35, graphinit=true, fg=fg)
# tree, _, _ = solveTree!(fg, tree);

# fg = RoME.generateCanonicalFG_Honeycomb!(42, graphinit=true, fg=fg)
# tree, _, _ = solveTree!(fg, tree);


## test numerical results

# working before IIF 1010

@warn("Test for beehive graph is using loose bounds until IIF #1010 is resolved.")
@test isapprox( getPPE(fg, :l11).suggested , [5;10*sin(pi/3)], atol=6)
@test isapprox( getPPE(fg, :l0).suggested , [20;0], atol=2)
@test isapprox( getPPE(fg, :l7).suggested , [20;-20*sin(pi/3)], atol=6)

# likely to fail until IIF 1010 is completed
@test_broken  isapprox( getPPE(fg, :x21).suggested[1:2] , [10;-20*sin(pi/3)], atol=4)


##

end

##

# using RoMEPlotting
# Gadfly.set_default_plot_size(40cm,25cm)

# ##

# plotSLAM2D(fg, drawContour=false, drawPoints=false, drawEllipse=false)

#