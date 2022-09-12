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
  nothing
end

