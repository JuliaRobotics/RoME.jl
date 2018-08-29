using RoME, IncrementalInference, TransformUtils, Distributions
using Base.Test
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

rr = RotationTest(MvNormal(zeros(3), 0.001*eye(3)))



@testset "Increased dimension root finding test" begin

# known rotations
eul = zeros(3,1)

R1 = zeros(3,1)
R2 = zeros(3,1)

res = randn(3)

rr(res, nothing, 1, (zeros(3,1),), R1, R2)
@test norm(res) < 1e-10



# random rotations
eul = 0.25*randn(3, 1)

R1 = rand(3,1)
R2 = rand(3,1)

res = zeros(3)

rr(res, nothing, 1, (zeros(3),), R1, R2)

@test norm(res) > 1e-3


end


@testset "test CommonConvWrapper functions" begin


N = 10

for i in 1:5


eul = 0.25*randn(3, N)
# res = zeros(3)
# @show rotationresidual!(res, eul, (zeros(0),x0))
# @show res
# gg = (res, x) -> rotationresidual!(res, eul, (zeros(0),x))
x0 = 0.1*randn(3)
res = zeros(3)
# @show gg(res, x0)
# @show res

A = rand(3,N)
B = rand(3,N)
At = deepcopy(A)
t = Array{Array{Float64,2},1}()
push!(t,A)
push!(t,B)
rr = RotationTest(MvNormal(zeros(3), 0.001*eye(3)))

zDim = 3

ccw = CommonConvWrapper(rr, t[1], zDim, t, measurement=(eul,))

@test ccw.xDim == 3

ccw(res, x0)

# and return complete fr/gwp
for n in 1:N

ccw.cpt[Threads.threadid()].particleidx = n
numericRootGenericRandomizedFnc!( ccw )

# test the result
qq = convert(Quaternion, Euler(eul[:,ccw.cpt[Threads.threadid()].particleidx]...))
q1 = convert(Quaternion, so3(ccw.cpt[Threads.threadid()].Y))
q2 = convert(Quaternion, so3(B[:,ccw.cpt[Threads.threadid()].particleidx]))
@test TransformUtils.compare(q1*q_conj(q2), qq, tol=1e-8)

end # particle for

end # i for

end # testset
