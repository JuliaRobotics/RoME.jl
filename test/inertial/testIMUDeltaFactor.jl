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

end

