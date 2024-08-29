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

include("variables/VariableTypes.jl")

end # module RoMETypes
