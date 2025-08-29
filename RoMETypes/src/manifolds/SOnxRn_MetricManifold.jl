struct MetricLieGroup{
    𝔽,
    O <: AbstractGroupOperation,
    M <: ManifoldsBase.AbstractManifold{𝔽},
    L <: LieGroup{𝔽, O, M},
    G <: RiemannianMetric,
} <: AbstractLieGroup{𝔽, O, M}
    lie_group::L
    metric::G
end

ManifoldsBase.base_manifold(G::MetricLieGroup) = G.lie_group
function ManifoldsBase.submanifold_component(G::MetricLieGroup, args...)
    return submanifold_component(G.lie_group, args...)
end
function ManifoldsBase.submanifold_components(G::MetricLieGroup, args...)
    return submanifold_components(G.lie_group, args...)
end
LieGroups.LieAlgebra(G::MetricLieGroup) = LieAlgebra(base_manifold(G))
LieGroups.inv!(G::MetricLieGroup, args...) = inv!(base_manifold(G), args...)
LieGroups.inv(G::MetricLieGroup, args...) = inv(base_manifold(G), args...)
LieGroups.compose!(G::MetricLieGroup, args...) = compose!(base_manifold(G), args...)
LieGroups.compose(G::MetricLieGroup, args...) = compose(base_manifold(G), args...)
function LieGroups.identity_element(G::MetricLieGroup, args...)
    return identity_element(base_manifold(G), args...)
end
function LieGroups.identity_element!(G::MetricLieGroup, args...)
    return identity_element!(base_manifold(G), args...)
end

# Left Invariant Rigid Body Kinematics Metric CrokeKumar eq 61.
# A family of left invariant metrics:
# G = [αI 0; 0 βI]
# where α and β are arbitrary constants, satisfies all the equations (80). 
# This are the only left-invariant metrics which are compatible with the acceleration connection.
# Can we add α and β as parameters to the metric?
struct LeftInvariantKinematicMetric <: RiemannianMetric end

function SOnxRn_MetricManifold(n)
    return MetricLieGroup(
        SpecialEuclideanGroup(n; variant = :right),
        LeftInvariantKinematicMetric(),
    )
end

SOnxRn_MetricManifoldType =
    Union{typeof(SOnxRn_MetricManifold(2)), typeof(SOnxRn_MetricManifold(3))}

# geodesics for metric (61) are the same as geodesics on the product manifold SO(3)×IR3
function Manifolds.exp(M::SOnxRn_MetricManifoldType, X)
    G = base_manifold(M)
    ε = identity_element(M, typeof(X))
    return exp(base_manifold(G), ε, X)
end

function Manifolds.exp!(M::SOnxRn_MetricManifoldType, g, X)
    G = base_manifold(M)
    ε = identity_element(M, typeof(g))
    return exp!(base_manifold(G), g, ε, X)
end

function ManifoldsBase.log(M::SOnxRn_MetricManifoldType, p)
    G = base_manifold(M)
    # ε = identity_element(G, typeof(p))
    # X = log(base_manifold(G), ε, p)
    PG = ProductLieGroup(map(LieGroup, G.manifold.manifolds, G.op.operations)...)
    X = log(PG, p)
    return X
end
function Manifolds.log!(M::SOnxRn_MetricManifoldType, X, p)
    G = base_manifold(M)
    ε = identity_element(M, typeof(p))
    log!(base_manifold(G), X, ε, p)
    return X
end

function Manifolds.inner(M::SOnxRn_MetricManifoldType, p, X, Y)
    Xtr = submanifold_components(M, X)[1]
    XRo = submanifold_components(M, X)[2]
    Ytr = submanifold_components(M, Y)[1]
    YRo = submanifold_components(M, Y)[2]
    # Metric on Chirikjian, p35 W = diagm([1,1,2])
    return dot(Xtr, Ytr) + dot(XRo, YRo) / 2
end

function DFG.getPointIdentity(::typeof(SOnxRn_MetricManifold(2)))
    return ArrayPartition(SA[0; 0.0], SA[1 0; 0 1.0])
end
function DFG.getPointIdentity(::typeof(SOnxRn_MetricManifold(3)))
    return ArrayPartition(SA[0, 0, 0.0], SA[1 0 0; 0 1 0; 0 0 1.0])
end

# LieGroups.LieAlgebra(G::SOnxRn_MetricManifoldType) = LieAlgebra(base_manifold(G))
