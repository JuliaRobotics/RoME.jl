# FIXME, deprecate old, refactor to better standard

export
  IIF,  # Aliases for various packages
  KDE,
  TU,
  AMP,
  DFG

# from Manifolds
export SpecialEuclidean


export
  # some transform functions
  cart2pol,
  pol2cart,

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

  Point3,
  Point3Point3,
  PackedPoint3Point3,

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
  PackedPose3Pose3XYYaw,
  PackedPriorPose3ZRP,
  Pose3Pose3Rotation,
  # Various utilities
  passTypeThrough,
  buildGraphChain!,

  # SLAM specific functions
  SLAMWrapper


export Pose3Pose3, PackedPose3Pose3
export VelPose2VelPose2, PackedVelPose2VelPose2
export PriorPose3, PackedPriorPose3



export generateGraph_ZeroPose, generateGraph_Hexagonal
export generateGraph_TwoPoseOdo, generateGraph_Circle
export generateGraph_Beehive!, generateGraph_Honeycomb!
export generateGraph_Boxes2D!
export generateGraph_Helix2D!, generateGraph_Helix2DSlew!, generateGraph_Helix2DSpiral!

export warmUpSolverJIT


export
  # camera model -- TODO --separate out
  CameraIntrinsic,
  CameraExtrinsic,
  CameraModelFull,
  project!,
  project,
  # keep
  cameraResidual!,

  # Didson convenience function
  addLinearArrayConstraint,

  # FG Analysis tools
  rangeErrMaxPoint2,
  rangeCompAllPoses,
  rangeCompAllPoses,

  # older features
  measureMeanDist,
  predictBodyBR,
  calcPosePointBearingRange,
  odomKDE,
  initFactorGraph!,
  addOdoFG!,
  getLastPose,
  getLastPose2D,
  addposeFG!,
  newLandm!,
  addBRFG!,
  addMMBRFG!,
  addAutoLandmBR!,
  projNewLandm!,
  malahanobisBR,
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

  # Feature tracking code
  Feature,
  initTrackersFrom,
  propAllTrackers!,
  measUpdateTrackers!,
  assocMeasWFeats!,

  lsrBR,

  # jld required Features Type
  LaserFeatures



# from odometry utils
export getMeasurementParametric
export accumulateDiscreteLocalFrame!, duplicateToStandardFactorVariable, extractDeltaOdo, resetFactor!
export odomKDE
export assembleChordsDict

export   DynPose2, DynPose2VelocityPrior, PackedDynPose2VelocityPrior, DynPose2Pose2, PackedDynPose2Pose2

# InertialPose3
export
  PreintegralCompensationGradients,
  InertialPose3Container,
  oplus,
  ⊕,
  getSample,
  InertialPose3,
  PackedInertialPose3,
  PriorInertialPose3,
  PackedPriorInertialPose3,
  compare