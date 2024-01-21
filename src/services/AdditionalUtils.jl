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
  
  # repeat for Pose3
  fg = initfg()
  addVariable!(fg, :x0, Pose3)
  addVariable!(fg, :x1, Pose3)
  addFactor!(fg, [:x0], PriorPose3(MvNormal(zeros(6), diagm(ones(6)))))
  addFactor!(fg, [:x0;:x1], Pose3Pose3(MvNormal(zeros(6), diagm(ones(6)))))
  
  initAll!(fg)

  IIF.initParametricFrom!(fg, :default)
  IIF.solveGraphParametric!(fg)

  solveGraph!(fg; multithread=true)

  nothing
end

"""
    $SIGNATURES

Return a PosePose factor type based on input homography matrix.
"""
function makePosePoseFromHomography(
  pHq::AbstractMatrix;
  covar=diagm([1,1,1,0.1,0.1,0.1].^2)
)
  len = size(pHq,1)-1
  M = SpecialEuclidean(len)
  e0 = ArrayPartition(SVector(0,0,0.), SMatrix{len,len}(1,0,0,0,1,0,0,0,1.))
  pTq = ArrayPartition(SVector(pHq[1:len,end]...), SMatrix{len,len}(pHq[1:len,1:len]))
  e0_Cpq = vee(M, e0, log(M, e0, pTq))
  if len == 2
    # use a covariance of 3x3 for Pose2 case
    covar_ = [covar[1,1] covar[1,2] covar[1,end]; covar[2,1] covar[2,2] covar[2,end]; covar[end,1] covar[end,2] covar[end,end]]
    Pose2Pose2(MvNormal(e0_Cpq, covar_))
  else
    Pose3Pose3(MvNormal(e0_Cpq, covar))
  end
end