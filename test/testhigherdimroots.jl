
using LinearAlgebra
using RoME
# , IncrementalInference, TransformUtils, Distributions
using Test
import  IncrementalInference: getSample

##

mutable struct RotationTest <: IncrementalInference.AbstractRelativeRoots
  z::MvNormal
end

getSample(cfo::CalcFactor{<:RotationTest}, N::Int=1) = (reshape(rand(cfo.factor.z,N),3,:),)

# 3 dimensional line, z = [a b][x y]' + c
function (cfo::CalcFactor{<:RotationTest})( meas, 
                                            var1, 
                                            var2)
  #
  #FIXME JT - I'm createing new res for simplicity, it may not hold up well though
  res = Vector{eltype(var1)}(undef, 3)
  
  z = meas
  dq = convert(Quaternion, Euler(z...))
  @show var1
  @show s1 = TU.so3(var1)
  s2 = TU.so3(var2)
  q1 = convert(Quaternion, s1)
  q2 = convert(Quaternion, s2)
  q12 = q1*q_conj(q2)
  qq = dq*q_conj(q12)
  @show res
  TransformUtils.vee!(res, convert(TU.so3, qq))
  return res
end

rr = RotationTest(MvNormal(zeros(3), 0.001*diagm(ones(3))))

##

@testset "Increased dimension root finding test" begin

## known rotations

eul = zeros(3)

R1 = zeros(3)
R2 = zeros(3)

res = randn(3)

# should be Sphere3
res = testFactorResidualBinary( rr, 
                                ContinuousEuclid{3}, 
                                ContinuousEuclid{3}, 
                                R1, 
                                R2, 
                                (eul,))
# rr(res, nothing, 1, (zeros(3,1),), R1, R2)


@test norm(res) < 1e-10

##

# random rotations
eul = 0.25*randn(3)

R1 = rand(3)
R2 = rand(3)

# res = zeros(3)

res = testFactorResidualBinary( rr, 
                                ContinuousEuclid{3}, 
                                ContinuousEuclid{3}, 
                                R1, 
                                R2, 
                                (eul,))
# rr(res, nothing, 1, (zeros(3),), R1, R2)


@test norm(res) > 1e-3

##

end


#