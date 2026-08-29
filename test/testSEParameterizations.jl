using LieGroups
using DistributedFactorGraphs
using Distributions
using StaticArrays
using DistributedFactorGraphs: ArrayPartition

T = ArrayPartition{Float64, Tuple{Matrix{Float64}, Vector{Float64}}}
G = SpecialEuclideanGroup(2; variant=:right)
ε = identity_element(G, T)

# T-R Coordinates, Chirikjian, p.35 Jr, above 10.76
function Jr(::Val{:SE2_TR}, Xⁱ::AbstractVector)
    @assert length(Xⁱ) == 3 "Xⁱ must be a vector of length 3"
    # Xⁱ = [x1, x2, θ]
    θ = Xⁱ[3]
    return [ 
         cos(θ)  sin(θ) 0; 
        -sin(θ)  cos(θ) 0;
          0           0 1
    ]
end

# Exponential Coordinates, Chirikjian, p.36 Jr, after 10.76
function Jr(::Val{:SE2_EXP}, Xⁱ::AbstractVector)
    @assert length(Xⁱ) == 3 "Xⁱ must be a vector of length 3"
    v1 = Xⁱ[1]
    v2 = Xⁱ[2]
    α = Xⁱ[3]
    if α == 0 
        return [1 0 0; 
                0 1 0; 
                0 0 1.]
    end
    return [
            (sin(α))/α   (1 - cos(α))/α   (α*v1 - v2 + v2*cos(α) - v1*sin(α))/α^2;
        (cos(α) - 1)/α         sin(α)/α   (v1 + α*v2 - v1*cos(α) - v2*sin(α))/α^2;
                     0                0                                         1
    ]
end

function Ad(::typeof(SpecialEuclideanGroup(2; variant=:right)), p)
  t = p.x[1]
  R = p.x[2]
  vcat(
      hcat(R, -SA[0 -1; 1 0]*t),
      SA[0 0 1]
  )
end

Ad(G, p)


##
# T-R Coordinates
X₀aⁱ = [10, 1, pi/4]
# LieAlgebra(G) is a Fiber?
X₀a = hat(LieAlgebra(G), X₀aⁱ, T)
# base_manifold(G) -> use the default Riemannian metric
p = exp(base_manifold(G), ε, X₀a)
# Exponential Coordinates
X₀b = log(G, p)

# J = XbJp * pJXa
# J =  jacobian_of_log(SE2, p)_wrt_p * jacobain_of_exp(SO2xEuclid2, Xa)_wrt_Xa
# J = J(Log(p), p) * J(exp(Xa), Xa)
J = inv(Jr(Val(:SE2_EXP), vee(LieAlgebra(G), X₀b))) * Jr(Val(:SE2_TR), vee(LieAlgebra(G), X₀a))

Σα = [
    1   0    0; 
    0 0.1    0; 
    0   0 0.01
]

N = 500000
# generate noise as tangent T-R Coordinates (assuming observation was in T-R)
Xⁱs = rand(MvNormal(Σα), N)

# map from T-R to exponential tangent representation
function φ(X) # work on better name
    # use Riemannian Exponential map to get the points (Chirikjian, p34, (eq. 10.75))
    np = LieGroups.compose(G, p, exp(base_manifold(G), ε, X)) # noisy p
    # convert the points back using the group logarithm map 
    return log(G, p, np) # -> Y
end

# Convert noisy coordinates in T-R coordinates to Exponential Coordinates at expansion point p
# is this a push-forward? disagreeement from Mateusz/Ronnie...
Yⁱs = map(eachcol(Xⁱs)) do Xⁱ
    # vector in a tangent space represented from T-R coordinates, [hybrid tangent representation]
    X = hat(LieAlgebra(G), Xⁱ, T) 
    # map from T-R to exponential tangent representation
    Y = φ(X)
    # get coordinates in exponential representation - Chirikjian, p35, (eq. 10.76)
    Yⁱ = vee(LieAlgebra(G), Y)
    return Yⁱ
end


Yⁱs  = hcat(Yⁱs...)

fit_Σβ = cov(fit_mle(MvNormal, Yⁱs))

# is fit_Σβ now the covariance converted from Σα in T-R coordinates to exponential coordinates at point p on G?

# is this a valid test?

# J = XbJp * pJXa
# J =  jacobian_of_log(SE2, p)_wrt_p * jacobain_of_exp(SO2xEuclid2, Xa)_wrt_Xa
# J = J(Log(p), p) * J(exp(Xa), Xa)
J = inv(Jr(Val(:SE2_EXP), vee(LieAlgebra(G), Y₀))) * Jr(Val(:SE2_TR), vee(LieAlgebra(G), X₀));

# Sola 2018, eq. 55
prop_Σβ = J * Σα * J'
fit_Σβ
# why do this? comparing what with what?
# What: linear propagate vs nonparametric mapping test of Covariance matrix
# Why: Trying resolve how to do robotics rigid transforms with new LieGroups.jl
isapprox(prop_Σβ, fit_Σβ; atol=1e-3)


## =================================================================================
##
## =================================================================================

# set Chirikjian 10.75 = 10.76
# x1 = (v2*(-1 + cos(α)) + v1*sin(α))/α
# x2 = (v1*( 1 - cos(α)) + v2*sin(α))/α
# θ = α

function map_TR_to_exp_coords(g)
    x1 = g[1]
    x2 = g[2]
    θ = g[3]

    α = θ
    
    #FIXME is this correct?
    if α == 0
        return [x1, x2, θ]
    end
    
    #FIXME maybe simplify this equations further
    v2 = (x1*α - x1*cos(α)*α - x2*sin(α)*α) / (2*(cos(α) - 1))
    v1 = (x2*α -  v2*sin(α)) / (1 - cos(α))

    return [v1, v2, α]
end

function map_exp_to_TR_coords(g)
    v1 = g[1]
    v2 = g[2]
    α = g[3]

    θ = α

    x1 = (v2*(-1 + cos(α)) + v1*sin(α))/α
    x2 = (v1*( 1 - cos(α)) + v2*sin(α))/α
  
    return [x1, x2, θ]
end

X₀ⁱ = [10, 1, pi/2]
X₀ = hat(LieAlgebra(G), X₀ⁱ, T)
p = exp(base_manifold(G), ε, X₀)
Y₀ = log(G, p)
Y₀ⁱ = vee(LieAlgebra(G), Y₀)

Yⁱ = map_TR_to_exp_coords(X₀ⁱ)
isapprox(Y₀ⁱ, Yⁱ)

Xⁱ = map_exp_to_TR_coords(Yⁱ)
isapprox(X₀ⁱ, Xⁱ)

N = 50000
# generate noise as tangent T-R Coordinates in the lie algebra
Xⁱs = rand(MvNormal(Σα), N)
Yⁱs = map(eachcol(Xⁱs)) do Xⁱ
    Yⁱ = map_TR_to_exp_coords(Xⁱ)
    return Yⁱ
end

Yⁱs  = hcat(Yⁱs...)

fit_Σβ = cov(fit_mle(MvNormal, Yⁱs))

Jfd = FiniteDiff.finite_difference_jacobian(map_TR_to_exp_coords, [0.,0,0])
prop_Σβ = Jfd * Σα * Jfd'
fit_Σβ

    X = hat(LieAlgebra(G), Xⁱ, T)
    q = exp(G, p, X)
    vee(LieAlgebra(G), log(G, p, q))
end

jac_Xbs = hcat(jac_Xbs...)


scatter(Xbs[1:2, :]; axis=(aspect=DataAspect(),), markersize=3)
scatter!(jac_Xbs[1:2, :]; markersize=3)