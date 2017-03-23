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
  HDF5,
  ProgressMeter


import Base: +
import TransformUtils: ⊖, ⊕, convert, compare
import IncrementalInference: convert, getSample, reshapeVec2Mat  #, compare


export
  # pass throughs from TransformUtils
  SE2,
  se2vee,
  se2vee!,
  SE3,
  Euler,
  Quaternion,
  AngleAxis,
  SO3,
  so3,
  compare,
  convert,


  # pass throughs from IncrementalInference
  ls,
  addFactor!,
  addNode!,
  getVert,
  getVertKDE,
  getVal,
  getData,
  FNDencode,
  FNDdecode,
  VNDencoder,
  VNDdecoder,
  GenericWrapParam,
  wipeBuildNewTree!,
  inferOverTree!,
  inferOverTreeR!,
  writeGraphPdf,
  savejld,
  loadjld,
  # overloaded functions from IIF
  # decodefg,
  # convertfrompackedfunctionnode,

  # RoME specific functions
  initfg,
  measureMeanDist,
  predictBodyBR,
  getLastPose,
  getLastPose2D,
  odomKDE,
  initFactorGraph!,
  addOdoFG!,
  addposeFG!,
  newLandm!,
  addBRFG!,
  addMMBRFG!,
  addAutoLandmBR!,
  projNewLandm!,
  malahanobisBR,

  # types
  BetweenPoses,

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
  +,
  evalPotential,
  getSample!,
  getSample,
  # obsolete
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

  # Point2D with null hypotheses
  PriorPoint2DensityNH,
  PackedPriorPoint2DensityNH,

  # acoustics
  Pose2DPoint2DBearingRangeDensity,
  PackedPose2DPoint2DBearingRangeDensity,
  Pose2DPoint2DRangeDensity,
  PackedPose2DPoint2DRangeDensity,

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
  PackedPriorPose3,
  Pose3Pose3,
  PackedPose3Pose3,
  projectParticles,
  ⊕,
  Pose3Pose3NH,
  PackedPose3Pose3NH,

  # partial Pose3
  PartialPriorRollPitchZ,
  PackedPartialPriorRollPitchZ,
  PartialPose3XYYaw,
  PackedPartialPose3XYYaw,

  # MultipleFeatures2D constraint functions
  MultipleFeatures2D,
  getUvecScaleFeature2D,
  getUvecScaleBaseline2D,

  # SLAM specific functions
  SLAMWrapper,

  # FG Analysis tools
  rangeErrMaxPoint2,
  rangeCompAllPoses,
  rangeCompAllPoses,


  # new robot navigation functionality
  triggerPose,
  GenericInSituSystem,
  InSituSystem,
  makeInSituSys,
  makeGenericInSituSys,
  advOdoByRules,
  progressExamplePlot,
  plotTrckStep,
  poseTrigAndAdd!,
  poseTrigAndAdd!,
  processTreeTrackersUpdates!


abstract BetweenPoses <: IncrementalInference.FunctorPairwise


include("BayesTracker.jl")

include("SensorModels.jl")
include("CameraModel.jl")
include("Point2D.jl")
include("Pose2D.jl")
include("BearingRange2D.jl")

include("BearingRangeDensity2D.jl")

include("Pose3Pose3.jl")
include("PartialPose3.jl")
include("MultipleFeaturesConstraint.jl")

include("Slam.jl")

include("RobotViz.jl")
include("RobotUtils.jl")

include("SimulationUtils.jl")

include("FactorGraphAnalysisTools.jl")

include("WheeledRobotUtils.jl")
include("NavigationSystem.jl")

include("fgos.jl")

# include("dev/ISAMRemoteSolve.jl")



end
