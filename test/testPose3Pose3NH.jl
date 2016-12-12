using Base.Test


begin
  function ff(X::Array{Float64,2})
    return X.+1.0
  end

  c = Categorical([0.5;0.5])
  A = zeros(3,5)

  skipdos = rand(c, 5)
  dos = skipdos.==2

  B = zeros(3,5)
  B[:,dos] = ff(A[:,dos])

  @test sum(sum(B[:,dos],2) .== sum(dos),1)[1,1] == 3
  println("Null Hypothesis matrix substitutions syntax works.")
end


warn("incomplete Pose3Pose3NH tests.")
