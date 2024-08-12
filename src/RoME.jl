module RoME

using Reexport

@reexport using IncrementalInference
@reexport using TransformUtils
@reexport using Distributions
@reexport using KernelDensityEstimate
@reexport using ApproxManifoldProducts

using
  Dates,
  FileIO,
  Distributed,
  LinearAlgebra,
  Statistics,
  CoordinateTransformations,
  # JLD2,
  ProgressMeter,
  DocStringExtensions,
  DistributedFactorGraphs,
  TensorCast,
  ManifoldsBase

using StaticArrays
using PrecompileTools

# to avoid name conflicts
import Manifolds
using Manifolds: hat, ProductGroup, ProductManifold, SpecialEuclidean, SpecialOrthogonal, TranslationGroup, identity_element, submanifold_component, Identity, affine_matrix

import Manifolds: project, project!, identity_element

import Rotations as _Rot
import Rotations: ⊕, ⊖ # TODO deprecate

export SpecialOrthogonal, SpecialEuclidean
export submanifold_component
# using Graphs,  # TODO determine how many parts still require Graphs still directly


import Base: +, \, convert
import TransformUtils: ⊖, ⊕, convert, ominus, veeQuaternion
import IncrementalInference: MB
import IncrementalInference: convert, getSample, reshapeVec2Mat, DFG
import IncrementalInference: getMeasurementParametric
import IncrementalInference: preambleCache
import IncrementalInference: InstanceType
# not sure why this is gives import error
import DistributedFactorGraphs: compare, @defVariable
import DistributedFactorGraphs: getDimension, getManifold

using OrderedCollections: OrderedDict
# const AMP = ApproxManifoldProducts



# export the API
include("ExportAPI.jl")


# load the source files
include("entities/SpecialDefinitions.jl")

#uses DFG v0.10.2 @defVariable for above
include("services/FixmeManifolds.jl")
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
include("factors/VelPosRotVelPos.jl")
include("factors/DynPose2D.jl")
include("factors/VelPose2D.jl")
include("factors/Point3D.jl")
include("factors/Point3Point3.jl")
include("factors/Pose3D.jl")
include("factors/Pose3Pose3.jl")
include("factors/PartialPose3.jl")
include("factors/MultipleFeaturesConstraint.jl")
include("factors/InertialPose3.jl")
# needs maintenance
include("factors/RangeAzimuthElevation.jl")
include("factors/Inertial/IMUDeltaFactor.jl")

# additional tools
include("services/FactorGraphAnalysisTools.jl")

# tools related to robotics
include("legacy/BayesTracker.jl")
include("factors/SensorModels.jl")
include("legacy/CameraModel.jl")
include("legacy/Slam.jl")
include("services/RobotUtils.jl")
include("services/ManifoldUtils.jl")
include("services/BearingRangeUtils.jl")
include("services/SimulationUtils.jl")
include("services/OdometryUtils.jl")
include("entities/RobotDataTypes.jl")
include("legacy/NavigationSystem.jl")

# generate canonical graphs
include("canonical/GenerateCommon.jl")
include("canonical/GenerateCircular.jl")
include("canonical/GenerateBox.jl")
include("canonical/GenerateHexagonal.jl")
include("canonical/GenerateHoneycomb.jl")
include("canonical/GenerateBeehive.jl")
include("canonical/GenerateHelix.jl")

# more utils requiring earlier functions
include("services/AdditionalUtils.jl")
include("services/g2oParser.jl")

# ScalarFields
include("services/ScalarFields.jl")

include("../ext/WeakdepsPrototypes.jl")
include("../ext/factors/GenericProjection.jl")
include("../ext/factors/InertialDynamic.jl")
include("../ext/factors/MixtureFluxPose2Pose2.jl")


# things on their way out
include("Deprecated.jl")


# manifold conversions required during transformation

@compile_workload begin
  # In here put "toy workloads" that exercise the code you want to precompile
  warmUpSolverJIT()
end

end
