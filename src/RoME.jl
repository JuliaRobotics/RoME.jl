module RoME

using
  IncrementalInference,
  Graphs,
  KernelDensityEstimate,
  Colors,
  Gadfly

export
  initfg,
  measureMeanDist,
  predictBodyBR,
  getLastPose2D,
  odomKDE,
  initFactorGraph!,
  addOdoFG!,
  newLandm!,
  addBRFG!,
  addMMBRFG!,
  addAutoLandmBR!,
  projNewLandm!,
  malahanobisBR,
  Feature,

  # helper functions
  getLastLandm2D,
  getLastPose2D,
  getNextLbl,

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
