
using Distributions, LinearAlgebra
using RoME
using Test


##

@testset "Test VelPos3 solve of of one variable" begin
##

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :w_P_1, RoME.VelPos3)
addFactor!(
  fg, [:w_P_1;], 
  PriorVelPos3(
    MvNormal(
      vcat(zeros(3), [100;200;300.]),
      diagm(ones(6))
    )
  )
  # ManifoldPrior(
  #   getManifold(RoME.VelPos3),
  #   ArrayPartition(zeros(3), [100;200;300.]),
  #   MvNormal(diagm(ones(6)))
  # )
)


stopping_criterion = 
  RoME.IncrementalInference.Manopt.StopAfterIteration(100) | 
  RoME.IncrementalInference.Manopt.StopWhenGradientNormLess(1e-12) | 
  RoME.IncrementalInference.Manopt.StopWhenStepsizeLess(1e-12)
debug = [:Iteration, :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :last_step_successful, "\n", :Stop]

IIF.solveGraphParametric!(
  fg; 
  stopping_criterion, 
  debug, 
  is_sparse=true, 
  damping_term_min=1e-12, 
  expect_zero_residual=true
)


## and serde

filepath = joinpath(tempdir(),"testgraphvepos")
saveDFG(filepath, fg)
fg_ = loadDFG(filepath)

@test length(ls(fg)) == length(ls(fg_))
@test length(lsf(fg)) == length(lsf(fg_))

##
end