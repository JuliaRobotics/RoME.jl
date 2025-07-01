

# Left Invariant Rigid Body Kinematics Metric CrokeKumar eq 61.
# A family of left invariant metrics:
# G = [αI 0; 0 βI]
# where α and β are arbitrary constants, satisfies all the equations (80). 
# This are the only left-invariant metrics which are compatible with the acceleration connection.
# Can we add α and β as parameters to the metric?
struct LeftInvariantKinematicMetric <: RiemannianMetric end

SOnxRn_MetricManifold(n) = MetricManifold(SpecialEuclideanGroup(n; variant=:right), LeftInvariantKinematicMetric())

SOnxRn_MetricManifoldType = Union{typeof(SOnxRn_MetricManifold(2)), typeof(SOnxRn_MetricManifold(3))}

# geodesics for metric (61) are the same as geodesics on the product manifold SO(3)×IR3
function Manifolds.exp(M::SOnxRn_MetricManifoldType, p, X)
    G = base_manifold(M)
    ε = identity_element(M) #, typeof(p))
    return LieGroups.compose(G, p, exp(base_manifold(G), ε, X))
end

function Manifolds.exp!(M::SOnxRn_MetricManifoldType, q, p, X)
    G = base_manifold(M)
    ε = identity_element(M) #, typeof(p))
    return LieGroups.compose!(G, q, p, exp(base_manifold(G), ε, X))
end

function Manifolds.log(M::SOnxRn_MetricManifoldType, p, q)
    G = base_manifold(M)
    ε = identity_element(M) #, typeof(p))
    X = log(base_manifold(G), ε, LieGroups.compose(G, inv(G, p), q))
    return X
end

function Manifolds.log!(M::SOnxRn_MetricManifoldType, X, p, q)
    G = base_manifold(M)
    ε = identity_element(M) #, typeof(p))
    log!(base_manifold(G), X, ε, LieGroups.compose(G, inv(G, p), q))
    return X
end

function Manifolds.inner(M::SOnxRn_MetricManifoldType, p, X, Y)
    Xtr = submanifold_components(M, X)[1]
    XRo = submanifold_components(M, X)[2]
    Ytr = submanifold_components(M, Y)[1]
    YRo = submanifold_components(M, Y)[2]
    # Metric on Chirikjian, p35 W = diagm([1,1,2])
    return dot(Xtr,Ytr) + dot(XRo, YRo)/2
end

function Manifolds.hat(M::SOnxRn_MetricManifoldType, p, X::AbstractVector, T=typeof(p))  
   return hat(LieAlgebra(base_manifold(M)), X, T)
end

function ManifoldsBase.get_vector(
    M::SOnxRn_MetricManifoldType, g, c, B::AbstractBasis{<:Any,TangentSpaceType}
)
    G = base_manifold(M)
    return get_vector(
        LieAlgebra(G),
        c,
        B;
        tangent_vector_type=ManifoldsBase.tangent_vector_type(G, typeof(g)),
    )
end

function ManifoldsBase.get_vector!(
    M::SOnxRn_MetricManifoldType, X, g, c, B::AbstractBasis{<:Any,TangentSpaceType}
)
    G = base_manifold(M)
    return get_vector!(LieAlgebra(G), X, c, B)
end

function ManifoldsBase.get_coordinates(
    M::SOnxRn_MetricManifoldType, g, X, B::AbstractBasis{<:Any,TangentSpaceType}
)
    G = base_manifold(M)
    return get_coordinates(LieAlgebra(G), X, B)
end

#FIXME: maybe replace this with identity_element
# call into base_manifold's identity_element, but got ambiguities.
Manifolds.identity_element(::typeof(SOnxRn_MetricManifold(2))) = ArrayPartition(SA[0;0.0],SA[1 0; 0 1.0])
Manifolds.identity_element(::typeof(SOnxRn_MetricManifold(3))) = ArrayPartition(SA[0,0,0.0],SA[1 0 0; 0 1 0; 0 0 1.0])

#FIXME remove, only temporary workaround as first signiture is used in many places
@deprecate LieGroups.identity_element(G::LieGroup, p)  identity_element(G, typeof(p)) false

# FIXME why is this still needed, hopefully can be removed soon 🐛💥
Base.convert(::Type{<:Tuple}, ::typeof(SOnxRn_MetricManifold(2))) = (:Euclid,:Euclid,:Circular)


## =======================================================================================
## SE2 + metric
## TODO move to test file
if false
M = SOnxRn_MetricManifold(2)
G = base_manifold(M)
ε = identity_element(M)
T = typeof(identity_element(M)) 
Xⁱ = [10, 1, pi/4]
X = hat(LieAlgebra(G), Xⁱ, T)
p = exp(base_manifold(G), ε, X)
q = compose(G, p, exp(base_manifold(G), ε, X))

q ≈ exp(M, p, X)
X ≈ log(M, p, q)


Xⁱ = [1, 1, 1]
X = hat(LieAlgebra(G), Xⁱ, T)
p = exp(base_manifold(G), ε, X)

W = diagm([1,1,2])

XX = hat(LieAlgebra(G), Xⁱ)
0.5*tr(XX*W*XX')
inner(M, p, X, X)
inner(base_manifold(G), p, X, X)

Manifolds.distance(M, ε, p)
Manifolds.distance(base_manifold(G), ε, p)

## SE3 + metric

M = SOnxRn_MetricManifold(3)
G = base_manifold(M)
ε = identity_element(M)
T = typeof(identity_element(M)) 
Xⁱ = [10, 1, 0, 0, 0, pi/4]
X = hat(LieAlgebra(G), Xⁱ, T)
p = exp(base_manifold(G), ε, X)
q = compose(G, p, exp(base_manifold(G), ε, X))

q ≈ exp(M, p, X)
X ≈ log(M, p, q)


Xⁱ = [0, 0, 0, 0, 1, 1]
X = hat(LieAlgebra(G), Xⁱ, T)
p = exp(base_manifold(G), ε, X)

0.5*tr(X.x[2]*X.x[2]')
dot(X.x[2], X.x[2])

inner(M, p, X, X)
inner(base_manifold(G), p, X, X)

Manifolds.distance(M, ε, p)
Manifolds.distance(base_manifold(G), ε, p)

Xⁱ = [1, 1, 0, 0, 1, 1]
W = diagm([1,1,1,2])
X = hat(LieAlgebra(G), Xⁱ, T)
XX = hat(LieAlgebra(G), Xⁱ)
0.5*tr(XX*W*XX')
inner(M, p, X, X)
inner(base_manifold(G), p, X, X)


end