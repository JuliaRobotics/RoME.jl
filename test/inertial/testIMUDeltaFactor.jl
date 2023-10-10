using Test
using Manifolds
using StaticArrays
using Rotations
using LinearAlgebra
using DistributedFactorGraphs
using RoME
using RoME: IMUDeltaGroup
using Dates

##
M = SpecialOrthogonal(3)
ΔR = exp(M, Identity(M), hat(M, Identity(M), [0.1, 0.2, 0.3]))
a = RotationVec(ΔR)
b = Rotations.AngleAxis(ΔR)


@testset "IMUDeltaFactor spot checks" begin

M = IMUDeltaGroup()
ϵ = identity_element(M)
@test ϵ == getPointIdentity(M)
@test inv(M, ϵ) == ϵ
@test compose(M, ϵ, ϵ) == ϵ

# vΔt, aΔt, ωΔt, Δt
X = hat(M, SA[0.0,0,0, 0,0,0, 0,0,1, 1] * 0.001)
p = exp(M, ϵ, X)
@test log(M, p) ≈ X

A = zeros(5,5)
A[3:4,5] .= 0.001 
exp(A)

affine_matrix(M, p)

Xc = SA[0.01, 0.02, 0.03,   0, 0, 0,   0.1, 0.2, 0.3,   1] * 0.001
X = hat(M, Xc)
p = exp(M, ϵ, X)
@test log(M, p) ≈ X
@test vee(M, log(M, p)) ≈ Xc

p = ArrayPartition(SMatrix{3,3}(1.0I), SA[1.,0,0], SA[0.,0,0], 0.0)
q = ArrayPartition(SMatrix{3,3}(1.0I), SA[1.,0,0], SA[0.1,0,0], 0.1)
# z-up (gravity positive) so free falling in negative z direction.
Δpq = RoME.boxminus(M, p, q;  g⃗ = SA[0,0,9.81])
@test Δpq.x[1] ≈ SMatrix{3,3}(1.0I)
@test Δpq.x[2] ≈ [0, 0, -9.81*0.1]
@test Δpq.x[3] ≈ [0, 0, -0.5*9.81*0.1^2] 
@test Δpq.x[4] == 0.1

# z-down (gravity negative) so free falling in positive z direction.
Δpq = RoME.boxminus(M, p, q;  g⃗ = SA[0,0,-9.81])
@test Δpq.x[1] ≈ SMatrix{3,3}(1.0I)
@test Δpq.x[2] ≈ [0, 0, 9.81*0.1]
@test Δpq.x[3] ≈ [0, 0, 0.5*9.81*0.1^2] 
@test Δpq.x[4] == 0.1


# vΔt, aΔt, ωΔt, Δt
X = hat(M, SA[0,0,0, 0,0,1.0, 0,0,0.5, 1] * 0.01)
p = exp(M, ϵ, X)
@test isapprox(p, ArrayPartition([1 -0.005 0.0; 0.005 1 0.0; 0 0 1], [0, 0, 0.01], [0, 0, 5.0e-5], 0.01), atol=1e-4)

θ=asin(0.1)*10 # for precicely 0.1
X = hat(M, SA[1,0,0, 0,0,0, 0,0,θ, 1] * 0.1)
p = exp(M, ϵ, X)

M_SE3 = SpecialEuclidean(3)
X_SE3 = hat(M_SE3, getPointIdentity(M_SE3), SA[1,0,0, 0,0,θ] * 0.1)
p_SE3 = exp_lie(M_SE3, X_SE3)
@test isapprox(p.x[3], p_SE3.x[1])

## test factor with rotation around z axis and initial velocity up
dt = 0.1
σ_a = 1e-4#0.16e-3*9.81  # noise density m/s²/√Hz
σ_ω = deg2rad(0.0001)  # noise density rad/√Hz
gn = MvNormal(diagm(ones(3)*σ_ω^2 * 1/dt))
an = MvNormal(diagm(ones(3)*σ_a^2 * 1/dt))

Σy  = diagm([ones(3)*σ_a^2; ones(3)*σ_ω^2])
gyros = [SA[0, 0, 0.001] + rand(gn) for _ = 1:11]
accels = [SA[0, 0, 9.81 - 1] + rand(an) for _ = 1:11]
timestamps = collect(range(0; step=dt, length=11))

a_b = SA[0.,0,0]
ω_b = SA[0.,0,0]

fac = RoME.IMUDeltaFactor(
    accels,
    gyros,
    timestamps,
    Σy,
    a_b,
    ω_b
)

# Rotation part
M_SO3 = SpecialOrthogonal(3)
ΔR = identity_element(M_SO3)
for g in gyros[1:end-1]
    exp!(M_SO3, ΔR, ΔR, hat(M_SO3, Identity(M_SO3), g*dt))
end
#TODO I would have expected these 2 to be exactly the same
@test isapprox(M_SO3, fac.Δ.x[1], ΔR; atol=1e-5)
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
Δpq = RoME.boxminus(RoME.IMUDeltaGroup(), p, q; g⃗ = SA[0,0,-9.81])
#TODO confirm gravity sign
@test Δpq.x[1] ≈ ΔR
@test Δpq.x[2] ≈ [0, 0, 8.81]
@test Δpq.x[3] ≈ [0, 0, 8.81/2] 
@test Δpq.x[4] == 1.0

end

