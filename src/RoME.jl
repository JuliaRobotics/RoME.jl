module RoME

using Reexport

@reexport using IncrementalInference
@reexport using TransformUtils
@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts

using
  Distributed,
  LinearAlgebra,
  Statistics,
  Graphs,  # TODO determine how many parts still require Graphs still directly
  Rotations,
  CoordinateTransformations,
  JLD2,
  ProgressMeter,
  DocStringExtensions

import Base: +, \, convert
import TransformUtils: ⊖, ⊕, convert, compare, ominus, veeQuaternion
import IncrementalInference: convert, getSample, reshapeVec2Mat, extractdistribution, DFG

# const AMP = ApproxManifoldProducts

export
  IIF,  # Aliases for various packages
  KDE,
  TU,
  AMP,
  # initfg,
  # RoME specific functions
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
  # getAll2D,
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
  PackedPoint2DPoint2D,
  Point2Point2WorldBearing,
  PackedPoint2Point2WorldBearing,
  Point2Point2Range,
  PackedPoint2Point2Range,
  PriorPoint2,
  PackedPriorPoint2,
  Pose2Point2BearingRange,
  Pose2Point2BearingRangeMH,
  PackedPose2Point2BearingRange,
  PackedPose2Point2BearingRangeMH,
  Pose2Point2Bearing,
  PackedPose2Point2Bearing,
  Pose2Point2Range,
  PriorPoint2,
  PackedPriorPoint2,
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

  # Polar types
  Polar,
  PolarPolar,
  PriorPolar,

  # MultipleFeatures2D constraint functions
  MultipleFeatures2D,
  getUvecScaleFeature2D,
  getUvecScaleBaseline2D,

  # Pose3, Three dimensional
  Pose3,
  Point3,
  PriorPoint3,
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
  basicFactorGraphExample,

  # jld required Features Type
  LaserFeatures,

  # Deprecated
  PartialPriorRollPitchZ,
  PackedPartialPriorRollPitchZ,
  PartialPose3XYYaw,
  PackedPartialPose3XYYaw,
  PartialPose3XYYawNH,
  PackedPartialPose3XYYawNH
  # Point2DPoint2DRange,
  # PackedPoint2DPoint2DRange,
  # PackedPoint2DPoint2DRange,
  # PackedPriorPoint2D,
  # Pose2DPoint2DRangeDensity, # to be deprecated
  # Point2DPoint2D, # deprecated
  # Point2DPoint2DRange, # deprecated
  # PriorPoint2D, # deprecated
  # Pose2DPoint2DBearingRange, # begin deprecated
  # Pose2DPoint2DBearing, # deprecated
  # PriorPoint2D, # deprecated
  # PackedPriorPoint2D # deprecated`


# doesnt seem to work
# @info "Setting IncrementalInference de-serialization namespace RoME"
# setSerializationNamespace!("RoME" => RoME)
# @info "done..."


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

# 2D
include("variables/Point2D.jl")
include("variables/Pose2D.jl")
include("variables/DynPoint2D.jl")
include("variables/DynPose2D.jl")

# 3D
include("variables/Point3D.jl")
include("variables/Pose3D.jl")


include("factors/Point2D.jl")
include("factors/Polar.jl")
include("factors/Pose2D.jl")
include("factors/Bearing2D.jl")
include("factors/Range2D.jl")
include("factors/BearingRange2D.jl")
include("factors/DynPoint2D.jl")
include("factors/DynPose2D.jl")
include("factors/Point3D.jl")
include("factors/Pose3Pose3.jl")
include("factors/PartialPose3.jl")
include("factors/MultipleFeaturesConstraint.jl")
include("factors/InertialPose3.jl")

include("Slam.jl")

include("RobotUtils.jl")

include("SimulationUtils.jl")

include("FactorGraphAnalysisTools.jl")

include("RobotDataTypes.jl") #WheeledRobotUtils
include("NavigationSystem.jl")


include("Deprecated.jl")

# include("dev/ISAMRemoteSolve.jl")



end
