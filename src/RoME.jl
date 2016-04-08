module RoME

using
  IncrementalInference,
  KernelDensityEstimate

export
  initfg,
  measureMeanDist,
  predictBodyBR,
  getLastPose2D,
  odomKDE,
  addOdoFG!,
  newLandm!,
  addBRFG!,
  addMMBRFG!,
  projNewLandm!,
  malahanobisBR,
  initFactorGraph!,

  # Some vizualization tools
  plotLsrScanFeats,
  drawFeatTrackers,
  saveImgSeq,
  lsrBR



include("BayesTracker.jl")
include("RobotViz.jl")
include("RobotUtils.jl")
include("SimulationUtils.jl")

end
