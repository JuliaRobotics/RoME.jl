
using Test
using RoME



@testset "test endless cycle case, issue #754" begin


fg = generateCanonicalFG_Circle(8;graphinit=false)
                                  
getSolverParams(fg).graphinit = false
getSolverParams(fg).treeinit = true
getSolverParams(fg).limititers = 100

# smtasks = Task[]
# tree, smt, hist = solveTree!(fg; smtasks=smtasks, verbose=true, timeout=50, recordcliqs=ls(fg));


getSolverParams(fg).drawtree = true
getSolverParams(fg).showtree = true
# getSolverParams(fg).dbg = true


mkpath(getLogPath(fg))
fid = open(joinLogPath(fg,"csmVerbose.log"), "w")
smtasks = Task[]
tree, smt, hist = solveTree!(fg; smtasks=smtasks, verbose=true, verbosefid=fid, timeout=100, recordcliqs=ls(fg), injectDelayBefore=[3=>(dwnInitSiblingWaitOrder_StateMachine=>5)]);

# after finished
close(fid)
open(joinLogPath(fg, "csmLogicalReconstructMax.log"),"w") do io
  IIF.reconstructCSMHistoryLogical(getLogPath(fg), fid=io)
end

hists = fetchCliqHistoryAll!(smtasks);

open(joinLogPath(fg, "csmLogicalSequential.log"),"w") do io
  printCliqHistorySequential(hists, nothing, io)
end

open(joinLogPath(fg, "csmLogicalLogical.log"),"w") do io
  printCSMHistoryLogical(hists, io)
end



end




#