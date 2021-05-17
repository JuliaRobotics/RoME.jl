# beehive larger test

using Test
using RoME
using Distributed

while nprocs() < 4
  addprocs(1)
end

using RoME
@everywhere using RoME

##

@testset "sequential beehive, repeated growing and solving" begin

##

fg = RoME.generateCanonicalFG_Beehive!(7, graphinit=true)
tree, _, _ = solveTree!(fg);

fg = RoME.generateCanonicalFG_Beehive!(14, graphinit=true, fg=fg)
tree, _, _ = solveTree!(fg, tree);

tree_ = deepcopy(tree); fg_ = deepcopy(fg);
fg = RoME.generateCanonicalFG_Beehive!(21, graphinit=true, fg=fg)
tree, _, _ = solveTree!(fg  , tree);


# fg = RoME.generateCanonicalFG_Beehive!(28, graphinit=true, fg=fg)
# tree, _, _ = solveTree!(fg, tree);

# fg = RoME.generateCanonicalFG_Beehive!(35, graphinit=true, fg=fg)
# tree, _, _ = solveTree!(fg, tree);

# fg = RoME.generateCanonicalFG_Beehive!(42, graphinit=true, fg=fg)
# tree, _, _ = solveTree!(fg, tree);


## test numerical results

# working before IIF 1010

# @test isappox( getPPE(fg, :l11).suggested , [], atol= )

# likely to fail until IIF 1010 is completed


##


end

##

using RoMEPlotting
Gadfly.set_default_plot_size(40cm,25cm)

##

plotSLAM2D(fg, drawContour=false, drawPoints=false, drawEllipse=false)

#