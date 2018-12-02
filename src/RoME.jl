module RoME

using Reexport

@reexport using IncrementalInference
@reexport using TransformUtils
@reexport using Distributions
@reexport using KernelDensityEstimate

using
  Distributed,
  LinearAlgebra,
  Statistics,
  Graphs,
  Rotations,
  CoordinateTransformations,
  JLD2,
  ProgressMeter,
  DocStringExtensions

import Base: +, \, convert
import TransformUtils: ⊖, ⊕, convert, compare, ominus, veeQuaternion
import IncrementalInference: convert, getSample, reshapeVec2Mat, extractdistribution  #, compare


export
  initfg,  # RoME specific functions
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
  Point2,
  Point2Point2,
  Point2DPoint2D, # deprecated
  PackedPoint2DPoint2D,
  Point2Point2WorldBearing,
  PackedPoint2Point2WorldBearing,
  Point2Point2Range,
  PackedPoint2Point2Range,
  Point2DPoint2DRange, # deprecated
  PackedPoint2DPoint2DRange,
  PriorPoint2,
  PackedPriorPoint2,
  PriorPoint2D, # deprecated
  PackedPriorPoint2D,
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
  PriorPoint2D, # deprecated
  PackedPriorPoint2,
  PackedPriorPoint2D, # deprecated`
  # Point2D with null hypotheses
  PriorPoint2DensityNH,
  PackedPriorPoint2DensityNH,

  # Velocity in Point2 types
  DynPoint2,
  DynPoint2VelocityPrior,
  DynPoint2DynPoint2,
  VelPoint2VelPoint2,
  Point2Point2Velocity,
  PackedDynPoint2VelocityPrior,
  PackedVelPoint2VelPoint2,

  # likely to be deprecated
  solveLandm,
  solvePose2,
  solveSetSeps,
  addPose2Pose2!,


  # acoustics
  Pose2Point2BearingRangeDensity,
  PackedPose2Point2BearingRangeDensity,
  Pose2Point2RangeDensity,
  Pose2DPoint2DRangeDensity, # to be deprecated
  PackedPose2Point2RangeDensity,

  # Pose2D
  Pose2,
  PriorPose2,
  PackedPriorPose2,
  PartialPriorYawPose2,
  PackedPartialPriorYawPose2,
  Pose2Pose2,
  PackedPose2Pose2,
  # velocity in Pose2
  DynPose2,
  DynPose2VelocityPrior,
  PackedDynPose2VelocityPrior,
  VelPose2VelPose2,
  PackedVelPose2VelPose2,
  DynPose2Pose2,
  PackedDynPose2Pose2,
  # Will be deprecated
  addPose2Pose2,

  # MultipleFeatures2D constraint functions
  MultipleFeatures2D,
  getUvecScaleFeature2D,
  getUvecScaleBaseline2D,

  # Pose3, Three dimensional
  Pose3,
  Point3,
  # Prior, # moved to IIF
  PriorPose3,
  PackedPriorPose3,
  Pose3Pose3,
  PackedPose3Pose3,
  projectParticles,
  ⊕,
  Pose3Pose3NH,
  PackedPose3Pose3NH,

  # partial Pose3
  PriorPose3ZRP,
  Pose3Pose3XYYaw,
  # Various utilities
  passTypeThrough,

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
  LaserFeatures,

  IIF,
  KDE,
  TU,

  # Deprecated
  PartialPriorRollPitchZ,
  PackedPartialPriorRollPitchZ,
  PartialPose3XYYaw,
  PackedPartialPose3XYYaw,
  PartialPose3XYYawNH,
  PackedPartialPose3XYYawNH


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
include("DynPose2D.jl")
include("Pose3D.jl")
include("factors/BearingRange2D.jl")

# include("BearingRangeDensity2D.jl")

include("factors/Pose3Pose3.jl")
include("factors/PartialPose3.jl")
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
