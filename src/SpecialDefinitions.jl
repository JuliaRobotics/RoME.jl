
const IIF = IncrementalInference
const KDE = KernelDensityEstimate

@compat const VoidUnion{T} = Union{Void, T}

@compat const CTs = CoordinateTransformations
@compat const TUs = TransformUtils


# TODO remove in RoME v0.1.7+, for use with IIF v0.3.9 and beyond
# import IncrementalInference: SamplableBelief
# SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity, IIF.AliasingScalarSampler}
SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity}

# @compat abstract type BetweenPoses <: IncrementalInference.FunctorPairwise end


vectoarr2(v) = reshape(v, length(v),1)
