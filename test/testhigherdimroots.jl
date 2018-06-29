using RoME, IncrementalInference, TransformUtils, Distributions
using Base.Test


println("Increased dimension root finding test")

mutable struct RotationTest <: IncrementalInference.FunctorPairwise
  z::MvNormal
end

# 3 dimensional line, z = [a b][x y]' + c
function (rt::RotationTest)(res::Vector{Float64}, userdata, idx, meas, var1,var2)
  z = meas[1]
  q1 = convert(Quaternion, Euler(z...))
  s = so3(var2[:,idx])
  q2 = convert(Quaternion, s)
  qq = q1*q_conj(q2)
  vee!(res, convert(so3, qq))
  nothing
end


N = 1

for i in 1:10
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

  gwp = GenericWrapParam{RotationTest}(rr, t, 1, 1)

  @time gwp(res, x0)

  # gwp.activehypo
  # gwp.hypotheses
  # gwp.params

  @show gwp.varidx
  gwp.measurement = (eul, )
  @show zDim = 3
  fr = FastRootGenericWrapParam{RotationTest}(gwp.params[gwp.varidx], zDim, gwp)

  @test fr.xDim == 3

  # and return complete fr/gwp
  @time for gwp.particleidx in 1:N
    # gwp(x, res)
    numericRootGenericRandomizedFnc!( fr )

    # test the result
    @show q1 = convert(Quaternion, Euler(eul[gwp.particleidx]...))
    @show q2 = convert(Quaternion, so3(fr.Y))
    @test TransformUtils.compare(q1, q2, tol=1e-8)
  end

  # y = numericRootGenericRandomizedFnc(
  #         gg,
  #         3, 3, x0   )

end
