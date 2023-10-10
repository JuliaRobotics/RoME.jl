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
Δpq = RoME.boxminus(M, p, q)
#TODO confirm test values
@test Δpq.x[1] ≈ SMatrix{3,3}(1.0I)
@test Δpq.x[2] ≈ [0, 0, -9.81*0.1]
@test Δpq.x[3] ≈ [0, 0, -0.5*9.81*0.1^2] 
@test Δpq.x[4] == 0.1

end

