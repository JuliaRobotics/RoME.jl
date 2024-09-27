using Test
using Manifolds
using StaticArrays
using Rotations
using LinearAlgebra
using DistributedFactorGraphs
using RoME
using RoME: IMUDeltaGroup
using Dates
using StaticArrays
# using ManifoldDiff

##
M = SpecialOrthogonal(3)
ΔR = exp(M, Identity(M), hat(M, Identity(M), [0.1, 0.2, 0.3]))
a = RotationVec(ΔR)
b = Rotations.AngleAxis(ΔR)

##
@testset "IMUDeltaFactor spot checks" begin
##

M = IMUDeltaGroup()
ϵ = identity_element(M)
@test ϵ == getPointIdentity(M)
@test inv(M, ϵ) == ϵ
@test compose(M, ϵ, ϵ) == ϵ

# vΔt, aΔt, ωΔt, Δt
X = hat(M, SA[0.0,0,0, 0,0,0, 0,0,1, 1] * 0.001)
p = exp(M, ϵ, X)
@test log(M, p) ≈ X

Xc = SA[0.1, 0.2, 0.3,  0.4, 0.5, 0.6,  0.7, 0.8, 0.9,  1] * 0.1
X = hat(M, Xc)
p = affine_matrix(M, exp(M, X))
affine_matrix(M, log(M, exp(M, X)))

X_af = RoME.vector_affine_matrix(M, X)
p_af = exp(X_af)
log(p_af)

@test isapprox(p, p_af)

# Xc = SVector{10,Float64}(vcat([0., 0, 0], aδt, ωδt, δt))

Xc = SA[0.01, 0.02, 0.03,   0, 0, 0,   0.1, 0.2, 0.3,   1] * 0.001
X = hat(M, Xc)
p = exp(M, ϵ, X)
@test log(M, p) ≈ X
@test vee(M, log(M, p)) ≈ Xc

Xc = SA[0, 0, 0,  0.01, 0.02, 0.03,  0.1, 0.2, 0.3,   1] * 0.001
X = hat(M, Xc)
p = exp(M, ϵ, X)
@test log(M, p) ≈ X
@test vee(M, log(M, p)) ≈ Xc

p = ArrayPartition(SMatrix{3,3}(1.0I), SA[1.,0,0], SA[0.,0,0], 0.0)
q = ArrayPartition(SMatrix{3,3}(1.0I), SA[1.,0,0], SA[0.1,0,0], 0.1)
# z-up (gravity positive)
Δpq = RoME.boxminus(M, p, q;  g⃗ = SA[0,0,9.81])
@test Δpq.x[1] ≈ SMatrix{3,3}(1.0I)
@test Δpq.x[2] ≈ [0, 0, 9.81*0.1]
@test Δpq.x[3] ≈ [0, 0, 0.5*9.81*0.1^2] 
@test Δpq.x[4] == 0.1

# z-down (gravity negative)
Δpq = RoME.boxminus(M, p, q;  g⃗ = SA[0,0,-9.81])
@test Δpq.x[1] ≈ SMatrix{3,3}(1.0I)
@test Δpq.x[2] ≈ [0, 0, -9.81*0.1]
@test Δpq.x[3] ≈ [0, 0, -0.5*9.81*0.1^2] 
@test Δpq.x[4] == 0.1


# vΔt, aΔt, ωΔt, Δt
X = hat(M, SA[0,0,0, 0,0,1.0, 0,0,0.5, 1] * 0.01)
p = exp(M, ϵ, X)
@test isapprox(p, ArrayPartition([1 -0.005 0.0; 0.005 1 0.0; 0 0 1], [0, 0, 0.01], [0, 0, 5.0e-5], 0.01), atol=1e-4)
X_af = RoME.vector_affine_matrix(M, X)
p_af = exp(X_af)
@test isapprox(affine_matrix(M, p), p_af, atol=1e-4)


# vΔt, aΔt, ωΔt, Δt
X = hat(M, SA[0,0,0, 1,0,0.0, 0,0,0, 1] * 0.01)
p = exp(M, ϵ, X)
@test isapprox(p, ArrayPartition([1.0 0 0; 0 1 0; 0 0 1], [0.01, 0, 0], [5e-5, 0, 0], 0.01), atol=1e-4)
X_af = RoME.vector_affine_matrix(M, X)
p_af = exp(X_af)
@test isapprox(affine_matrix(M, p), p_af, atol=1e-4)


X = hat(M, SA[1,0,0, 1,0,0.0, 0,0,0.01, 1])
p = exp(M, ϵ, X)
q = compose(M, p, exp(M, ϵ, X))
isapprox(compose(M, p, exp(M, ϵ, X)), exp(M, p, X))

RoME.adjointMatrix(M, X) * vee(M,X)

X_af = RoME.vector_affine_matrix(M, X)
p_af = affine_matrix(M, p)

Y = p_af*X_af*inv(p_af)
vee(M, ArrayPartition(Y[1:3,1:3], Y[1:3,4], Y[1:3,5], Y[4,5]))

#testing adjoint matrix with properties
Adₚ = RoME.AdjointMatrix(M, p)

q1 = compose(M, p, exp(M, X))
q2 = compose(M, exp(M, hat(M, Adₚ*vee(M, X))), p)
@test isapprox(q1, q2)

@test isapprox(RoME.AdjointMatrix(M, inv(M, p)), inv(Adₚ))

@test isapprox(
    RoME.vector_affine_matrix(M, hat(M, Adₚ*vee(M, X))),
    affine_matrix(M, p) * RoME.vector_affine_matrix(M, X) * affine_matrix(M, inv(M, p))
)

ad = RoME.adjointMatrix(M, X)
@test isapprox(exp(ad), Adₚ)


X = hat(M, SA[0.1, 0.2, 0.3,  0.4, 0.5, 0.6,  0.7, 0.8, 0.9,  1] * 0.1)
ad = RoME.adjointMatrix(M, X)
Adₚ = RoME.AdjointMatrix(M, exp(M, X))
@test isapprox(exp(ad), Adₚ)

Y = hat(M, SA[0.9, 0.8, 0.7,  0.6, 0.5, 0.4,  0.3, 0.2, 0.1,  1] * 0.1)
#TODO test Jr with something like
# Z1 = ManifoldDiff.differential_exp_argument_lie_approx(M, p, X, Y)
Z2 = RoME.Jr(M, X) * vee(M, Y)

#test right jacobian with [Ad(g)] = Jl*Jr⁻¹ - Chirikjian p29
jr = RoME.Jr(M, X)
jl = RoME.Jr(M, -X)
@test isapprox(jl*inv(jr), Adₚ)

θ=asin(0.1)*10 # for precisely 0.1
X = hat(M, SA[1,0,0, 0,0,0, 0,0,θ, 1] * 0.1)
p = exp(M, ϵ, X)

M_SE3 = SpecialEuclidean(3; vectors=HybridTangentRepresentation())
X_SE3 = hat(M_SE3, getPointIdentity(M_SE3), SA[1,0,0, 0,0,θ] * 0.1)
p_SE3 = exp_lie(M_SE3, X_SE3)
@test isapprox(p.x[3], p_SE3.x[1])

## test factor with rotation around z axis and initial velocity up
# DUPLICATED IN testInertialDynamic.jl
dt = 0.1
N = 11

σ_a = 1e-4 #0.16e-3*9.81  # noise density m/s²/√Hz
σ_ω = deg2rad(0.0001)  # noise density rad/√Hz
imu = RoME.generateField_InertialMeasurement_noise(; dt, N, rate=SA[0, 0, 0.001], accel0=SA[0, 0, 9.81-1], σ_a, σ_ω)

timestamps = collect(range(0; step=dt, length=N))

a_b = SA[0.,0,0]
ω_b = SA[0.,0,0]

fac = RoME.IMUDeltaFactor(
    imu.accels,
    imu.gyros,
    timestamps,
    imu.Σy,
    a_b,
    ω_b
)

# Rotation part
M_SO3 = SpecialOrthogonal(3)
ΔR = identity_element(M_SO3)
#NOTE internally integration is done from 2:end
# at t₀, we start at identity
# - the measurement at t₀ is from t₋₁ to t₀ and therefore part of the previous factor
# - the measurement at t₁ is from t₀ to t₁ with the time dt of t₁ - t₀
for g in imu.gyros[2:end]
    exp!(M_SO3, ΔR, ΔR, hat(M_SO3, Identity(M_SO3), g*dt))
end
@test isapprox(M_SO3, fac.Δ.x[1], ΔR)
#TODO The velocity and position parts might be off from expected because of the connection/metric of the group
# Velocity part
@test isapprox(fac.Δ.x[2], [0,0,8.81], atol=1e-3) # after 1 second at 9.81 m/s
# position part
@test isapprox(fac.Δ.x[3], [0,0,8.81/2], atol=1e-3) # after 1 second at 1/2*9.81*1^2 m
# Δt part
@test isapprox(fac.Δ.x[4], 1.0) 

#Same for delta variables 
p = ArrayPartition(SMatrix{3,3}(1.0I), SA[0.,0,0], SA[1.,0,0], 0.0)
q = ArrayPartition(SMatrix{3,3}(ΔR),   SA[0.,0,-1], SA[1.,0,-0.5], 1.0)
# z-down (gravity acc negative) so free falling in positive z direction.
Δpq = RoME.boxminus(RoME.IMUDeltaGroup(), p, q; g⃗ = SA[0,0,9.81])
#TODO confirm gravity sign
@test Δpq.x[1] ≈ ΔR
@test Δpq.x[2] ≈ [0, 0, 8.81]
@test Δpq.x[3] ≈ [0, 0, 8.81/2] 
@test Δpq.x[4] == 1.0

Σy  = diagm(ones(6)*0.1^2)
a_b = SA[0.,0,0]
ω_b = SA[0.,0,0]

fg = initfg()
fg.solverParams.graphinit = false

foreach(enumerate(Nanosecond.(timestamps[[1,end]] * 10^9))) do (i, nanosecondtime)
    addVariable!(fg, Symbol("x",i-1), RotVelPos; nanosecondtime)
end

addFactor!(
    fg,
    [:x0],
    ManifoldPrior(
        getManifold(RotVelPos),
        ArrayPartition(
            SA[ 1.0 0.0 0.0; 
                0.0 1.0 0.0; 
                0.0 0.0 1.0], 
            SA[10.0, 0.0, 0.0], 
            SA[0.0, 0.0, 0.0]
        ),
        MvNormal(diagm(ones(9)*1e-3))
    )
)

addFactor!(fg, [:x0, :x1], fac)

@time IIF.solveGraphParametric!(fg; is_sparse=false, damping_term_min=1e-12, expect_zero_residual=true);
# @time IIF.solveGraphParametric!(fg; stopping_criterion, debug, is_sparse=false, damping_term_min=1e-12, expect_zero_residual=true);

getVariableSolverData(fg, :x0, :parametric).val[1] ≈ ArrayPartition(SA[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], SA[10.0, 0.0, 0.0], SA[0.0, 0.0, 0.0])
x1 = getVariableSolverData(fg, :x1, :parametric).val[1]
@test isapprox(SpecialOrthogonal(3), x1.x[1], ΔR, atol=1e-5)
@test isapprox(x1.x[2], [10, 0, -1], atol=1e-3)
@test isapprox(x1.x[3], [10, 0, -0.5], atol=1e-3)


dt = 0.01
N = 11
dT = (N-1)*dt
imu = RoME.generateField_InertialMeasurement(;dt,N,accel0=SA[0, 0, 9.81],rate=SA[0, 0, 0.1])
# gyros = [SA[0, 0, 0.1] for _ = 1:N]
# accels = [SA[0, 0, 9.81] for _ = 1:N]
timestamps = collect(range(0; step=dt, length=N))

Δ, Σ, J_b = RoME.preintegrateIMU(imu.accels, imu.gyros, timestamps, Σy, a_b, ω_b)
Σ = Σ[SOneTo(9),SOneTo(9)]

@test Δ.x[1] ≈ RotZ(0.1*dT)
@test Δ.x[2] ≈ [0, 0, 9.81*dT]
@test Δ.x[3] ≈ [0, 0, 0.5*9.81*dT^2] 
@test Δ.x[4] == dT

##
# imu = RoME.generateField_InertialMeasurement(;dt,N,accel0=SA[0, 0, 9.81],rate=SA[0.01, 0, 0])
gyros = [SA[0.01, 0, 0] for _ = 1:N]
accels = [SA[0, 0, 9.81] for _ = 1:N]
timestamps = collect(range(0; step=dt, length=N))

Δ, Σ, J_b = RoME.preintegrateIMU(accels, gyros, timestamps, Σy, a_b, ω_b)

@test Δ.x[1] ≈ RotX(0.01*dT)
#just checking sign
@test Δ.x[2][2] < 0
@test Δ.x[3][2] < 0

##
# imu = RoME.generateField_InertialMeasurement(;dt,N,accel0=SA[0, 0, 9.81],rate=SA[0, 0.01, 0])
gyros = [SA[0, 0.01, 0] for _ = 1:N]
accels = [SA[0, 0, 9.81] for _ = 1:N]
timestamps = collect(range(0; step=dt, length=N))

Δ, Σ, J_b = RoME.preintegrateIMU(accels, gyros, timestamps, Σy, a_b, ω_b)

@test Δ.x[1] ≈ RotY(0.01*dT)
#just checking sign
@test Δ.x[2][1] > 0
@test Δ.x[3][1] > 0

##
end

##
