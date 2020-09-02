
# continuously drive in a circle to demonstrate fixed lag

using RoME


fg = initfg()
getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
getSolverParams(fg).dbg = true


defaultFixedLagOnTree!(fg, 5, limitfixeddown=true)


tree = BayesTree()

generateCanonicalFG_Circle(6, fg=fg, offsetPoses=0, stopEarly=3, cyclePoses=6)
tree, smt, hists = solveTree!(fg, tree, recordcliqs=ls(fg));
drawTree(tree,filepath=joinLogPath(fg,"tree3.pdf"))


generateCanonicalFG_Circle(6, fg=fg, offsetPoses=3, stopEarly=6, cyclePoses=6)
tree, smt, = solveTree!(fg, tree);
drawTree(tree,filepath=joinLogPath(fg,"tree6.pdf"))


generateCanonicalFG_Circle(12, fg=fg, offsetPoses=6, stopEarly=9, cyclePoses=6)
tree, smt, = solveTree!(fg, tree);
drawTree(tree,filepath=joinLogPath(fg,"tree9.pdf"))


generateCanonicalFG_Circle(12, fg=fg, offsetPoses=9, stopEarly=12, cyclePoses=6)
tree, smt, = solveTree!(fg, tree);
drawTree(tree,filepath=joinLogPath(fg,"tree12.pdf"))


generateCanonicalFG_Circle(18, fg=fg, offsetPoses=12, stopEarly=15, cyclePoses=6)
tree, smt, = solveTree!(fg, tree);
drawTree(tree,filepath=joinLogPath(fg,"tree15.pdf"))


generateCanonicalFG_Circle(18, fg=fg, offsetPoses=15, stopEarly=18, cyclePoses=6)
tree, smt, = solveTree!(fg, tree);
drawTree(tree,filepath=joinLogPath(fg,"tree18.pdf"))


drawGraph(fg, show=true)


## visualization of numerical results

using RoMEPlotting

Gadfly.set_default_plot_size(35cm,20cm)
plotSLAM2D(fg)
