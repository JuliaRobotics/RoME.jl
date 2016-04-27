module RoME

using
  IncrementalInference,
  Graphs,
  KernelDensityEstimate,
  Distributions,
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

  # helper functions
  get2DSamples,
  getAll2D,
  get2DSampleMeans,
  getAll2DMeans,
  getAll2DPoses,
  get2DPoseSamples,
  get2DPoseMeans,
  getKDE,
  getVertKDE,
  get2DPoseMax,
  getAll2DLandmarks,
  get2DLandmSamples,
  get2DLandmMeans,
  get2DLandmMax,

  # helper functions
  getLastLandm2D,
  getLastPose2D,
  getNextLbl,

  # some transform functions
  cart2pol,
  pol2cart,

  # Feature tracking code
  Feature,
  initTrackersFrom,
  propAllTrackers!,
  measUpdateTrackers!,
  assocMeasWFeats!,

  # Some vizualization tools
  plotLsrScanFeats,
  drawFeatTrackers,
  saveImgSeq,
  lsrBR,

  # draw pose beliefs etc
  drawPoses,
  drawLandms,
  drawPosesLandms,
  drawSubmaps,
  investigatePoseKDE




include("BayesTracker.jl")
include("RobotViz.jl")
include("RobotUtils.jl")
include("SimulationUtils.jl")

end
