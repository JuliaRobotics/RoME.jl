

# Rigid Body Kinematics Metric CrokeKumar eq 61.
# A family of left invariant metrics:
# G = [αI 0; 0 βI]
# where α and β are arbitrary constants, satisfies all the equations (80). 
# This are the only left-invariant metrics which are compatible with the acceleration connection.
# Can we add α and β as parameters to the metric?
struct LeftInvariantRBKMetric <: RiemannianMetric end

SOnxRnManifold(n) = MetricManifold(SpecialEuclideanGroup(n; variant=:right), LeftInvariantRBKMetric())

SOnxRnManifoldType = Union{typeof(SOnxRnManifold(2)), typeof(SOnxRnManifold(3))}

Manifolds.identity_element(::typeof(SOnxRnManifold(2))) = ArrayPartition(SA[0;0.0],SA[1 0; 0 1.0])
Manifolds.identity_element(::typeof(SOnxRnManifold(3))) = ArrayPartition(SA[0,0,0.0],SA[1 0 0; 0 1 0; 0 0 1.0])

# geodesics for metric (61) are the same as geodesics on the product manifold SO(3)×IR3
function Manifolds.exp(M::SOnxRnManifoldType, p, X)
    G = base_manifold(M)
    ε = identity_element(M)
    return compose(G, p, exp(base_manifold(G), ε, X))
end

function Manifolds.log(M::SOnxRnManifoldType, p, q)
    G = base_manifold(M)
    ε = identity_element(M)
    X = log(base_manifold(G), ε, compose(G, inv(G, p), q))
    return X
end

function Manifolds.inner(M::SOnxRnManifoldType, p, X, Y)
    Xtr = submanifold_components(M, X)[1]
    XRo = submanifold_components(M, X)[2]
    Ytr = submanifold_components(M, Y)[1]
    YRo = submanifold_components(M, Y)[2]
    # Metric on Chirikjian, p35 W = diagm([1,1,2])
    return dot(Xtr,Ytr) + dot(XRo, YRo)/2
end








## =======================================================================================
## SE2 + metric

M = SOnxRnManifold(2)
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

M = SOnxRnManifold(3)
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


