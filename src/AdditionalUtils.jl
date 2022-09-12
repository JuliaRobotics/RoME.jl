# Precompile Utilities

"""
    $SIGNATURES

Load and solve a canonical or user factor graph to warm up---precompile---several RoME/Caesar related functions.
"""
function warmUpSolverJIT(;fg::AbstractDFG=generateGraph_Hexagonal(),
                          drawtree::Bool=false )
  #
  
  fcts = ls(fg, :x0)
  fcts = ls(fg)
  fcts = lsf(fg, :x0f1)
  fcts = lsf(fg)
  drawGraph(fg, show=false)
  # getSolverParams(fg).drawtree = drawtree
  initAll!(fg)
  tree = solveTree!(fg)
  
  IIF.initParametricFrom!(fg, :default)
  IIF.solveGraphParametric!(fg)
  
  fg = generateGraph_Hexagonal(graphinit=false)
  initAll!(fg)
  
  fg = initfg()
  addVariable!(fg, :x0, Pose3)
  addVariable!(fg, :x1, Pose3)
  addFactor!(fg, [:x0], PriorPose3(MvNormal(zeros(6), diagm(ones(6)))))
  addFactor!(fg, [:x0;:x1], Pose3Pose3(MvNormal(zeros(6))))
  
  solveGraph!(fg)

  nothing
end

