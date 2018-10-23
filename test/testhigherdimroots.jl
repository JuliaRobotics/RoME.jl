
using LinearAlgebra
using RoME
# , IncrementalInference, TransformUtils, Distributions
using Test
# import  IncrementalInference: getSample


mutable struct RotationTest <: IncrementalInference.FunctorPairwise
  z::MvNormal
end

# 3 dimensional line, z = [a b][x y]' + c
function (rt::RotationTest)(res::Vector{Float64}, userdata, idx, meas, var1, var2)
  z = view(meas[1],:,idx)
  dq = convert(Quaternion, Euler(z...))
  s1 = so3(var1[:,idx])
  s2 = so3(var2[:,idx])
  q1 = convert(Quaternion, s1)
  q2 = convert(Quaternion, s2)
  q12 = q1*q_conj(q2)
  qq = dq*q_conj(q12)
  vee!(res, convert(so3, qq))
  nothing
end

global rr = RotationTest(MvNormal(zeros(3), 0.001*Matrix{Float64}(LinearAlgebra.I, 3,3)))



@testset "Increased dimension root finding test" begin

# known rotations

global eul = zeros(3,1)

global R1 = zeros(3,1)
global R2 = zeros(3,1)

global res = randn(3)

rr(res, nothing, 1, (zeros(3,1),), R1, R2)
@test norm(res) < 1e-10



# random rotations
global eul = 0.25*randn(3, 1)

global R1 = rand(3,1)
global R2 = rand(3,1)

global res = zeros(3)

rr(res, nothing, 1, (zeros(3),), R1, R2)

@test norm(res) > 1e-3


end


@testset "test CommonConvWrapper functions" begin


global N = 10

global i=1
for i in 1:5


global eul = 0.25*randn(3, N)

# res = zeros(3)
# @show rotationresidual!(res, eul, (zeros(0),x0))
# @show res
# gg = (res, x) -> rotationresidual!(res, eul, (zeros(0),x))
global x0 = 0.1*randn(3)
global res = zeros(3)
# @show gg(res, x0)
# @show res

global A = rand(3,N)
global B = rand(3,N)
global At = deepcopy(A)
global t = Array{Array{Float64,2},1}()
push!(t,A)
push!(t,B)
global rr = RotationTest(MvNormal(zeros(3), 0.001*Matrix{Float64}(LinearAlgebra.I, 3,3)))

global zDim = 3

global ccw = CommonConvWrapper(rr, t[1], zDim, t, measurement=(eul,)) # old bug where measurement is not patched through fixed in IIF v0.3.9

@test ccw.xDim == 3

# TODO remove
ccw.measurement = (eul,)

ccw(res, x0)

# and return complete fr/gwp
global n = 1
for n in 1:N

ccw.cpt[Threads.threadid()].particleidx = n
numericRootGenericRandomizedFnc!( ccw )

# test the result
global qq = convert(Quaternion, Euler(eul[:,ccw.cpt[Threads.threadid()].particleidx]...))
global q1 = convert(Quaternion, so3(ccw.cpt[Threads.threadid()].Y))
global q2 = convert(Quaternion, so3(B[:,ccw.cpt[Threads.threadid()].particleidx]))
@test TransformUtils.compare(q1*q_conj(q2), qq, tol=1e-8)

end # particle for

end # i for

end # testset
