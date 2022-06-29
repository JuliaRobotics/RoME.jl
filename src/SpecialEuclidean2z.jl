# for PoSpecialEuclidean2z which includes z height

import ManifoldsBase: check_point, check_vector, manifold_dimension, exp!, inner, representation_size, get_embedding
# import Manifolds: activate_traits
import Manifolds: submanifold_component, submanifold_components
import Manifolds: affine_matrix, _padpoint!, _padvector!
import Base: show

const SpecialEuclidean2z = SemidirectProductGroup{
    ℝ,
    TranslationGroup{Tuple{3},ℝ},
    SpecialOrthogonal{2},
    RotationAction{TranslationGroup{Tuple{3},ℝ},SpecialOrthogonal{2},LeftAction},
}

function SpecialEuclidean2z()
    Tz = TranslationGroup(3)
    SO2 = SpecialOrthogonal(2)
    A = RotationAction(Tz, SO2)
    return SemidirectProductGroup(Tz, SO2, A)
end

const SpecialEuclidean2zOperation = SemidirectProductOperation{
    RotationAction{TranslationGroup{Tuple{3},ℝ},SpecialOrthogonal{2},LeftAction},
}
const SpecialEuclidean2zIdentity = Identity{SpecialEuclidean2zOperation}

Base.show(io::IO, ::SpecialEuclidean2z) = print(io, "SpecialEuclidean2z()")

@inline function active_traits(f, M::SpecialEuclidean, args...)
    return merge_traits(IsGroupManifold(M.op), IsExplicitDecorator())
end

Base.@propagate_inbounds function submanifold_component(
    ::SpecialEuclidean2z,
    p::AbstractMatrix,
    ::Val{1},
) where {n}
    # think this is translation xyz from affine matrix?
    return view(p, 1:3, 2 + 1)
end

Base.@propagate_inbounds function submanifold_component(
    ::SpecialEuclidean2z,
    p::AbstractMatrix,
    ::Val{2},
) where {n}
    return view(p, 1:2, 1:2)
end

function submanifold_components(
    G::SpecialEuclidean2z,
    p::AbstractMatrix,
) where {n}
    @assert size(p) == (n + 1, n + 1)
    @inbounds t = submanifold_component(G, p, Val(1))
    @inbounds R = submanifold_component(G, p, Val(2))
    return (t, R)
end

Base.@propagate_inbounds function _padpoint!(
    ::SpecialEuclidean2z,
    q::AbstractMatrix,
)
    for i in 1:2
        q[3 + 1, i] = 0
    end
    q[3 + 1, 2 + 1] = 1
    return q
end

Base.@propagate_inbounds function _padvector!(
    ::SpecialEuclidean2z,
    X::AbstractMatrix,
)
    for i in 1:(2 + 1)
        X[3 + 1, i] = 0
    end
    return X
end

function affine_matrix(G::SpecialEuclidean2z, p)
    pis = submanifold_components(G, p)
    pmat = allocate_result(G, affine_matrix, pis...)
    map(copyto!, submanifold_components(G, pmat), pis)
    @inbounds _padpoint!(G, pmat)
    return pmat
end
affine_matrix(::SpecialEuclidean2z, p::AbstractMatrix) = p

function check_point(G::SpecialEuclidean2z, p::AbstractMatrix; kwargs...)
    errs = DomainError[]
    # homogeneous
    if !isapprox(p[end, :], [zeros(size(p, 2) - 1)..., 1]; kwargs...)
        push!(
            errs,
            DomainError(
                p[end, :],
                "The last row of $p is not homogeneous, i.e. of form [0,..,0,1].",
            ),
        )
    end
    # translate part
    err2 = check_point(submanifold(G, 1), p[1:3, end]; kwargs...)
    !isnothing(err2) && push!(errs, err2)
    # SOn
    err3 = check_point(submanifold(G, 2), p[1:2, 1:2]; kwargs...)
    !isnothing(err3) && push!(errs, err3)
    
    if length(errs) > 1
        return CompositeManifoldError(errs)
    end
    return length(errs) == 0 ? nothing : first(errs)
end

function exp_lie!(G::SpecialEuclidean{2}, q, X)
    SO2 = submanifold(G, 2)
    b, Ω = submanifold_components(G, X)
    t, R = submanifold_components(G, q)
    @assert size(R) == (2, 2)
    @assert size(t) == (2,)
    @assert size(b) == (2,)

    θ = vee(SO2, identity_element(SO2, R), Ω)[1]
    sinθ, cosθ = sincos(θ)
    if θ ≈ 0
        α = 1 - θ^2 / 6
        β = θ / 2
    else
        α = sinθ / θ
        β = (1 - cosθ) / θ
    end

    @inbounds begin
        R[1] = cosθ
        R[2] = sinθ
        R[3] = -sinθ
        R[4] = cosθ
        t[1] = α * b[1] - β * b[2]
        t[2] = α * b[2] + β * b[1]
        t[3] = b[3] + 
        _padpoint!(G, q)
    end
    return q
end