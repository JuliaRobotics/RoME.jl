module RoME

using Reexport

@reexport using IncrementalInference
@reexport using TransformUtils
@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts

using
  Dates,
  Distributed,
  LinearAlgebra,
  Statistics,
  Graphs,  # TODO determine how many parts still require Graphs still directly
  Rotations,
  CoordinateTransformations,
  JLD2,
  ProgressMeter,
  DocStringExtensions,
  DistributedFactorGraphs,
  TensorCast

import Base: +, \, convert
import TransformUtils: ⊖, ⊕, convert, compare, ominus, veeQuaternion
import IncrementalInference: convert, getSample, reshapeVec2Mat, extractdistribution, DFG, getManifolds
# not sure why this is gives import error
import DistributedFactorGraphs: compare
import DistributedFactorGraphs: getDimension, getManifolds
# const AMP = ApproxManifoldProducts

const InstanceType{T} = Union{Type{<:T},T}

export
  IIF,  # Aliases for various packages
  KDE,
  TU,
  AMP,
  DFG,

  # RoME specific functions
  measureMeanDist,
  predictBodyBR,
  calcPosePointBearingRange,
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
  get2DSampleMeans,
  getAll2DMeans,
  getAll2DPoses,
  get2DPoseSamples,
  get2DPoseMeans,
  get2DPoseMax,
  get2DLandmSamples,
  get2DLandmMeans,
  get2DLandmMax,

  # helper functions
  getLastLandm2D,
  getLastPose2D,
  getNextLbl,

  # RobotUtils
  getRangeKDEMax2D,
  getLastPoses,
  setSolvableOldPoses!,

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
  # backprojectRandomized!,
  residual!,
  residualLRBE!,
  reuseLBRA,
  ominus,
  ominus!,
  +,
  evalPotential,
  getSample!,
  getSample,

  # Didson convenience function
  addLinearArrayConstraint,

  # camera model -- TODO --separate out
  CameraIntrinsic,
  CameraExtrinsic,
  CameraModelFull,
  project!,
  project,
  # keep
  cameraResidual!,

  # Point2D
  Point2,
  Point2Point2,
  PackedPoint2Point2,
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
  PackedPose2Point2Range,
  PriorPoint2,
  PackedPriorPoint2,

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
  compareDensity,

  # Pose2D
  Pose2,
  PriorPose2,
  PackedPriorPose2,
  PartialPriorYawPose2,
  PackedPartialPriorYawPose2,
  Pose2Pose2,
  PackedPose2Pose2,
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
  predictVariableByFactor,

  # jld required Features Type
  LaserFeatures,

  # Deprecated
  PartialPriorRollPitchZ,
  PackedPartialPriorRollPitchZ,
  PartialPose3XYYaw,
  PackedPartialPose3XYYaw


include("SpecialDefinitions.jl")

## More variable types
# 2D
# include("variables/Point2D.jl")
# include("variables/Pose2D.jl")
# include("variables/DynPoint2D.jl")
# include("variables/DynPose2D.jl")

# 3D
# include("variables/Point3D.jl")
# include("variables/Pose3D.jl")

#uses DFG v0.10.2 @defVariable for above
include("variables/VariableTypes.jl")

## More factor types
# RoME internal factors (FYI outside factors are easy, see Caesar documentation)
include("factors/Point2D.jl")
include("factors/Range2D.jl")
include("factors/Bearing2D.jl")
include("factors/BearingRange2D.jl")
include("factors/Polar.jl")
include("factors/PriorPose2.jl")
include("factors/PartialPriorPose2.jl")
include("factors/Pose2D.jl")
include("factors/Pose2Point2.jl")
include("factors/MutablePose2Pose2.jl")
include("factors/DynPoint2D.jl")
include("factors/VelPoint2D.jl")
include("factors/DynPose2D.jl")
include("factors/VelPose2D.jl")
include("factors/Point3D.jl")
include("factors/Pose3D.jl")
include("factors/Pose3Pose3.jl")
include("factors/PartialPose3.jl")
include("factors/MultipleFeaturesConstraint.jl")
include("factors/InertialPose3.jl")

# tools that come and go
include("TemporaryFunctionality.jl")

# additional tools
include("FactorGraphAnalysisTools.jl")

# tools related to robotics
include("BayesTracker.jl")
include("SensorModels.jl")
include("CameraModel.jl")
include("Slam.jl")
include("RobotUtils.jl")
include("SimulationUtils.jl")
include("OdometryUtils.jl")
include("RobotDataTypes.jl")
include("NavigationSystem.jl")
include("CanonicalGraphs.jl")
include("g2oParser.jl")


# things on their way out
include("Deprecated.jl")

# optional tools
using Requires

function __init__()
  # combining neural networks natively into the non-Gaussian  factor graph object
  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin
    @info "RoME is adding Flux related functionality."
    include("factors/flux/models/Pose2OdoNN_01.jl") # until a better way is found to deserialize
    include("factors/flux/MixtureFluxPose2Pose2.jl")
    # include("factors/flux/FluxModelsPose2Pose2.jl")
  end
end

end
