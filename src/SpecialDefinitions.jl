
const IIF = IncrementalInference
const KDE = KernelDensityEstimate

const NothingUnion{T} = Union{Nothing, T}

const CTs = CoordinateTransformations
const TUs = TransformUtils



vectoarr2(v) = reshape(v, length(v),1)
