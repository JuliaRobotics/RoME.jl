module RoME

using
  IncrementalInference,
  Graphs,
  TransformUtils,
  KernelDensityEstimate,
  Distributions,
  Colors,
  Gadfly,
  JLD,
  HDF5

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

  # DIDSON sonar model
  LinearRangeBearingElevation,
  project!,
  project,
  backprojectRandomized!,
  residual!,
  ominus,
  evalPotential,
  getSample,

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
  drawMarginalContour,
  accumulateMarginalContours,

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
  investigatePoseKDE,

  # solve with isam in pytslam
  doISAMSolve,
  drawCompPosesLandm,

  # Victoria Park data specific
  LaserFeatures,
  addLandmarksFactoGraph!,
  appendFactorGraph!,
  doBatchRun,
  rotateFeatsToWorld,

  togglePrtStbLines




include("BayesTracker.jl")
include("RobotViz.jl")
include("RobotUtils.jl")
include("SensorModels.jl")
include("SimulationUtils.jl")
include("VictoriaParkTypes.jl")
include("VicPrkEstimator.jl")
include("dev/ISAMRemoteSolve.jl")

end
