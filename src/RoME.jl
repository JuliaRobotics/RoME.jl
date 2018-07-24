module RoME

using
  IncrementalInference,
  Graphs,
  TransformUtils,
  CoordinateTransformations,
  Rotations,
  KernelDensityEstimate,
  Distributions,
  JLD,
  HDF5,
  ProgressMeter,
  DocStringExtensions,
  Compat

import Base: +, \, convert
import TransformUtils: ⊖, ⊕, convert, compare, ominus, veeQuaternion
import IncrementalInference: convert, getSample, reshapeVec2Mat, extractdistribution  #, compare


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
  FunctorSingleton,
  FunctorPairwise,
  FunctorPairwiseNH,   # will become obsolete
  FunctorSingletonNH,  # will become obsolete
  ls,
  addFactor!,
  addNode!,
  getVert,
  getVertKDE,
  getVal,
  setVal!,
  getData,
  FNDencode,
  FNDdecode,
  localProduct,
  predictbelief,
  VNDencoder,
  VNDdecoder,
  GenericWrapParam,
  wipeBuildNewTree!,
  inferOverTree!,
  inferOverTreeR!,
  writeGraphPdf,
  savejld,
  loadjld,
  FactorGraph,
  initializeNode!,
  isInitialized,
  ensureAllInitialized!,
  getPoints,
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
  veePose3,
  veePose,
  \,
  RangeAzimuthElevation,

  # types
  # BetweenPoses,

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

  # RobotUtils
  getRangeKDEMax2D,

  # some transform functions
  cart2pol,
  pol2cart,

  # Feature tracking code
  Feature,
  initTrackersFrom,
  propAllTrackers!,
  measUpdateTrackers!,
  assocMeasWFeats!,

  lsrBR,

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
  Point2DPoint2D,
  PackedPoint2DPoint2D,
  Point2Point2WorldBearing,
  PackedPoint2Point2WorldBearing,
  Point2DPoint2DRange,
  PackedPoint2DPoint2DRange,
  PriorPoint2,
  PriorPoint2D, # to be deprecated
  PackedPriorPoint2,
  PackedPriorPoint2D, # to be deprecated`
  DynPoint2VelocityPrior,
  DynPoint2DynPoint2,
  VelPoint2VelPoint2,
  Point2Point2Velocity,
  Pose2Point2BearingRange,
  Pose2DPoint2DBearingRange, # begin deprecated
  Pose2Point2BearingRangeMH,
  PackedPose2Point2BearingRange,
  PackedPose2Point2BearingRangeMH,
  Pose2Point2Bearing,
  Pose2DPoint2DBearing, # deprecated
  Pose2Point2Range,
  Point2DPoint2DRange,
  PackedPoint2DPoint2DRange,
  PriorPoint2,
  PriorPoint2D, # to be deprecated
  PackedPriorPoint2,
  PackedPriorPoint2D, # to be deprecated`
  solveLandm,
  solvePose2,
  solveSetSeps,
  addPose2Pose2!,

  # Point2D with null hypotheses
  PriorPoint2DensityNH,
  PackedPriorPoint2DensityNH,

  # acoustics
  Pose2Point2BearingRangeDensity,
  PackedPose2Point2BearingRangeDensity,
  Pose2Point2RangeDensity,
  Pose2DPoint2DRangeDensity, # to be deprecated
  PackedPose2Point2RangeDensity,

  Pose2,
  Point2,
  DynPoint2,
  Pose3,
  Point3,
  Prior,

  # Pose2D
  PriorPose2,
  PackedPriorPose2,
  # Pose2Pose2_NEW,
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
  PartialPose3XYYawNH,
  PackedPartialPose3XYYawNH,

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
  poseTrigAndAdd!,
  poseTrigAndAdd!,
  processTreeTrackersUpdates!,
  addSoftEqualityPoint2D,
  vectoarr2,

  # jld required Features Type
  LaserFeatures

  # # solve with isam in pytslam
  # doISAMSolve,
  # drawCompPosesLandm,
  #
  # # Victoria Park data specific
  # addLandmarksFactoGraph!,
  # appendFactorGraph!,
  # doBatchRun,
  # rotateFeatsToWorld





include("SpecialDefinitions.jl")

include("BayesTracker.jl")

include("SensorModels.jl")
include("CameraModel.jl")
include("Point2D.jl")
include("DynPoint2D.jl")
include("Pose2D.jl")
include("Pose3D.jl")
include("BearingRange2D.jl")

# include("BearingRangeDensity2D.jl")

include("Pose3Pose3.jl")
include("PartialPose3.jl")
include("MultipleFeaturesConstraint.jl")

include("InertialPose3.jl")

include("Slam.jl")

include("RobotUtils.jl")

include("SimulationUtils.jl")

include("FactorGraphAnalysisTools.jl")

include("RobotDataTypes.jl") #WheeledRobotUtils
include("NavigationSystem.jl")


include("Deprecated.jl")

# include("dev/ISAMRemoteSolve.jl")



end
