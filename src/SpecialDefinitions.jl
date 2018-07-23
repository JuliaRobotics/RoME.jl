
SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity}

# @compat abstract type BetweenPoses <: IncrementalInference.FunctorPairwise end

@compat const VoidUnion{T} = Union{Void, T}

@compat const CTs = CoordinateTransformations
@compat const TUs = TransformUtils

vectoarr2(v) = reshape(v, length(v),1)
