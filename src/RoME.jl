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

import IncrementalInference: getSample

export
  initfg,
  measureMeanDist,
  predictBodyBR,
  getLastPose,
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

  togglePrtStbLines,

  # Didson model
  evalPotential,
  LinearRangeBearingElevation,
  project!,
  project,
  backprojectRandomized!,
  residual!,
  residualLRBE!,
  reuseLBRA,
  ominus,
  ominus!,
  evalPotential,
  getSample!,
  getSample,
  WrapParam,
  WrapParamArray,

  # Didson convenience function
  addLinearArrayConstraint,

  # camera model -- TODO --separate out
  CameraIntrinsic,
  CameraExtrinsic,
  CameraModelFull,
  project!,
  project,
  backprojectRandomized!,
  # keep
  cameraResidual!,

  # Point2D
  Pose2DPoint2DBearingRange,
  PackedPose2DPoint2DBearingRange,
  Pose2DPoint2DRange,
  Point2DPoint2DRange,
  PackedPoint2DPoint2DRange,
  PriorPoint2D,
  PackedPriorPoint2D,
  solveLandm,
  solvePose2,
  solveSetSeps,
  addPose2Pose2!,

  # Pose2D
  PriorPose2,
  PackedPriorPose2,
  Pose2Pose2,
  PackedPose2Pose2,
  addPose2Pose2,

  # packedTypes
  passTypeThrough,

  # Pose3D
  PriorPose3,
  Pose3Pose3,
  projectParticles,
  ⊕,
  Pose3Pose3NH,
  PackedPose3Pose3NH,

  # SLAM specific functions
  SLAMWrapper


  # Datalayer converters
  # VNDencoder,
  # VNDdecoder,
  # FNDencode,
  # FNDdecode





include("BayesTracker.jl")

include("SensorModels.jl")
include("CameraModel.jl")
include("Point2D.jl")
include("Pose2D.jl")
include("BearingRange2D.jl")
include("Pose3Pose3.jl")

include("Slam.jl")

include("RobotViz.jl")
include("RobotUtils.jl")

include("SimulationUtils.jl")


include("VictoriaParkTypes.jl")
include("VicPrkEstimator.jl")

include("dev/ISAMRemoteSolve.jl")



end
