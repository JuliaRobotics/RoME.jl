
const IIF = IncrementalInference
const KDE = KernelDensityEstimate

const VoidUnion{T} = Union{Void, T}

const CTs = CoordinateTransformations
const TUs = TransformUtils



vectoarr2(v) = reshape(v, length(v),1)
