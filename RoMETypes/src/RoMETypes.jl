module RoMETypes

using DistributedFactorGraphs
using DocStringExtensions
using Manifolds
using RecursiveArrayTools
using StaticArrays

import DistributedFactorGraphs: getVariableType, AbstractManifoldMinimize

export 
    Point2,
    Point3,
    Pose2,
    Pose3,
    Rotation3,
    RotVelPos,
    VelPos3,
    DynPoint2,
    DynPose2,
    projectCartesian

export 
    PriorPoint2,
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

include("variables/VariableTypes.jl")
include("factors/FactorTypes.jl")

end # module RoMETypes
