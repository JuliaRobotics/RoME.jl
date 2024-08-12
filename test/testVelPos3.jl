
using Distributions, LinearAlgebra
using RoME
using Test


##

@testset "Test VelPos3 solve of of one variable" begin
##

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :w_GNSS_1, RoME.VelPos3)
addFactor!(
  fg, [:w_GNSS_1;], 
  ManifoldPrior(
    getManifold(RoME.VelPos3),
    ArrayPartition(zeros(3), [100;200;300.]),
    MvNormal(diagm(ones(6)))
  )
)


stopping_criterion=StopAfterIteration(100) | StopWhenGradientNormLess(1e-12) | StopWhenStepsizeLess(1e-12)
debug = [:Iteration, :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :last_step_successful, "\n", :Stop]

IIF.solveGraphParametric!(fg; stopping_criterion,  debug, is_sparse=true, damping_term_min=1e-12, expect_zero_residual=true)


##
end