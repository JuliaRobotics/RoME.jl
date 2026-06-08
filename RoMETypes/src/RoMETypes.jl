module RoMETypes

using DistributedFactorGraphs
using DocStringExtensions
using LieGroups
using RecursiveArrayTools
using StaticArrays
using LinearAlgebra

using ManifoldsBase:
    submanifold_components, TangentSpaceType, AbstractBasis, RiemannianMetric
using Manifolds: MetricManifold

import ManifoldsBase
import Manifolds

export Point2,
    Point3,
    Pose2,
    Pose3,
    Rotation3,
    RotVelPos,
    VelPos3,
    DynPoint2,
    DynPose2

export PriorPoint2,
    PackedPriorPoint2,
    PriorPoint3,
    PackedPriorPoint3,
    Pose2Pose2,
    PackedPose2Pose2,
    Pose3Pose3,
    PackedPose3Pose3,
    PriorPose3,
    PackedPriorPose3,
    Point3Point3,
    PackedPoint3Point3,
    Point2Point2,
    PackedPoint2Point2,
    PriorPose2,
    PackedPriorPose2

export LeftInvariantMetricSE

include("manifolds/LeftInvariantMetricSE.jl")
include("variables/VariableTypes.jl")
include("factors/FactorTypes.jl")

end # module RoMETypes
