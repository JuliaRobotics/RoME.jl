
# ## RANDOM DEBUG DEV CODE BELOW
#
# hist = getCliqSolveHistory(tree, :x1)
#
# tree = solveTree!(fg, skipcliqids=[:x1;:x6;:x4;:x3], recordcliqs=[:x2;:x0])
#
# plotTreeUpMsgs(fg, tree, :x1, levels=1)
#
#
# cliq = getClique(tree, :x2)
#
# # OLD
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
# plotKDE(dwinmsgs[:x1][1], dims=[1;2], levels=2)
#
# # NEW
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
# plotKDE(dwinmsgs[:x1][1], dims=[1;2], levels=2)
#
#
#
# ## Check init message for x3
#
# cliq = getClique(tree, :x1)
#
# # OLD
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
# # NEW
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
#
#
# ## Check init message for x4
#
# cliq = getClique(tree, :x4)
#
# # OLD
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=false)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
# # NEW
# dwinmsgs = IIF.prepCliqInitMsgsDown!(fg, tree, getParent(tree,cliq)[1], cliq, dbgnew=true)
# plotKDE(dwinmsgs[:x3][1], dims=[1;2], levels=2)
#
#
#
#
# plotKDE(fg, :x1, dims=[1;2])
# plotPose(fg, :x1)
#
# plotKDE(fg, [:x0;:x1;:x2;:x4;:x6], dims=[1;2],levels=1)
#
#
# drawTree(tree, imgs=true)
#
#
#
# ## check contents
#
# cliq = getClique(tree, :x1)
#
# getCliqMsgsUp(cliq)
#
#
# #
#
# # DEBUG after solveCliq :x2
#
#
# # solve third clq
# smt, hist = solveCliq!(fg, tree, :x1, cliqHistories=hist, recordcliq=true)
#
#
# # Plot to see what is going on
# Gadfly.set_default_plot_size(35cm,20cm)
# plotKDE(fg, [:x0;:x1;:x2], dims=[1;2],levels=1)
#
# plotKDE(fg, :x1, dims=[1;2],levels=1)
#
#
#
# ## DEBUG where is down message :x3
#
# prnt = getClique(tree, :x3)
# getCliqInitUpMsgs(prnt)
#
#
# assignTreeHistory!(tree, hist)
# printCliqHistorySummary(tree, :x1)
#
#
# # think the issue is in here
# # 21:01:00.856  8   null         attemptCliqInitUp     false null | upsolved upsolved
#
#
#
# # somethings up with cliq x1
#
# csmc_8_test = getCliqSolveHistory(tree, :x1)[8][4]
# csmc_9_test = getCliqSolveHistory(tree, :x1)[9][4]
#
# # Just before the mistake
# drawTree(csmc_8_test.tree, show=true)
# prnt_8_test = getClique(csmc_8_test.tree, :x3)
# getCliqInitUpMsgs(prnt_8_test)
#
# # just after the mistake
# prnt_9_test = getClique(csmc_9_test.tree, :x3)
# getCliqInitUpMsgs(prnt_9_test)
#
#
# # so develop  based on sandbox step 8
# stuff = sandboxCliqResolveStep(tree,:x1,8)
# getCliqInitUpMsgs(getClique(stuff[4].tree,:x3))
#
#
#
# IIF.getCliqMsgsUp(tree,:x1)
#
#
# ## DEBUG
#
#
#
#
#
#
#
# ## DEBUG after solveCliq :x2 difference in down init cycle when :dontUseParentFactorsInitDown
#
# assignTreeHistory!(tree, hist)
# printCliqHistorySummary(tree, :x2)
#
# # WORKS
#
# delete!(getSolverParams(getCliqSolveHistory(tree, :x2)[15][4].dfg).devParams, :dontUseParentFactorsInitDown)
#
# good = sandboxCliqResolveStep(tree,:x2,15)
#
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x1, title="(X2,) sandbox 15 good,")
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x2, title="(X2,) sandbox 15 good,")
#
#
# # BREAKS
#
# getSolverParams(getCliqSolveHistory(tree, :x2)[15][4].dfg).devParams[:dontUseParentFactorsInitDown] = ""
#
# bad = sandboxCliqResolveStep(tree,:x2,15)
#
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x1, title="(X2,) sandbox 15 bad,")
# # plotPairVariables(csmcX2_ref[16][4].cliqSubFg, stuff[4].cliqSubFg, :x2, title="(X2,) sandbox 15 bad,")
#
# plotPairVariables(good[4].cliqSubFg, bad[4].cliqSubFg, :x2, title="(X2,) sandbox 15,")
# # plotPairPose2(good[4].cliqSubFg, bad[4].cliqSubFg, :x2, title="(X2,) sandbox 15,")
#
# Gadfly.set_default_plot_size(30cm,18cm)
#
# ## DEBUG
#
#
#
#
#
#
#
# #3 First difference, up msg from (:x2) cliq is 'wrong'
#
# Gadfly.set_default_plot_size(30cm,18cm)
# um1 = IIF.getCliqMsgsUp(tree, :x0)
# plotKDE(um1[:x1], dims=[1;2], levels=2)
#
#
#
# um2 = IIF.getCliqMsgsUp(tree, :x2)
# plotKDE(um2[:x1], dims=[1;2], levels=2)
#
#
#
#
# csmcX2_test = getCliqSolveHistory(tree, :x2)
#
# # csmcX2_test[15][4].cliqSubFg
#
# plotKDE(getKDE(csmcX2_test[17][4].cliqSubFg, :x1),dims=[1;2],levels=2,title="csmcX2_test[17]")
#
#
#
#
# ## SOMETHING changes between csmcX2_ref[16/17] and csmcX2_test[16/17]
#
# writeGraphPdf(csmcX2_ref[17][4].cliqSubFg, show=true)
# writeGraphPdf(csmcX2_test[17][4].cliqSubFg, show=true)
#
#
# plotKDE(csmcX2_ref[16][4].cliqSubFg, [:x1;:x2;:x3], dims=[1;2])
# plotKDE(csmcX2_test[16][4].cliqSubFg, [:x1;:x2;:x3], dims=[1;2])
#
#
#
# isInitialized(csmcX2_test[16][4].cliqSubFg, :x3)
#
#
#
#
#
#
# # before init cycle
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x1)
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x2)
# plotPairVariables(csmcX2_ref[16][4].cliqSubFg, csmcX2_test[16][4].cliqSubFg, :x3)
#
#
# # after solve
# plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x1)
# plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x2)
# plotPairVariables(csmcX2_ref[17][4].cliqSubFg, csmcX2_test[17][4].cliqSubFg, :x3)
