using RoME
using StaticArrays
using Test

##

@testset "Test (partial) Point2 priors on Pose2" begin
##

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :x0, Pose2)
addVariable!(fg, :x1, Pose2)
addVariable!(fg, :x2, Pose2)

odo2 = MvNormal(SA[1.0; 0.0; 0.0], diagm([1.0, 1.0, 0.1]).^2)
addFactor!(fg, [:x0, :x1], Pose2Pose2(odo2))
addFactor!(fg, [:x1, :x2], Pose2Pose2(odo2))

mv2 = MvNormal(SA[0.0; 0.0;], SA[1.0 0; 0 1])
f = addFactor!(fg, [:x0;], PriorPoint2(mv2))

approxConvBelief(fg, getLabel(f), :x0)

mv2 = MvNormal(SA[0.0; 2.0;], SA[1.0 0; 0 1])
f = addFactor!(fg, [:x2;], PriorPoint2(mv2))

IIF.solveGraphParametric!(fg)
p0 = getVal(fg, :x0, solveKey=:parametric)[1]
p1 = getVal(fg, :x1, solveKey=:parametric)[1]
p2 = getVal(fg, :x2, solveKey=:parametric)[1]

# driving east to west (along +y in world), therefore all headings = pi/2
@test isapprox(M, p0, ArrayPartition([0, 0.0], [0 -1.0; 1.0 0]), atol=1e-6)
@test isapprox(M, p1, ArrayPartition([0, 1.0], [0 -1.0; 1.0 0]), atol=1e-6)
@test isapprox(M, p2, ArrayPartition([0, 2.0], [0 -1.0; 1.0 0]), atol=1e-6)

# piggy back test on PPE for parametric
@test all( getPPE.(fg, [:x0;:x1;:x2], :parametric) .|> s->isapprox(s.suggested[3], pi/2; atol=1e-6) )

## test similar on nonparametric solve

solveGraph!(fg)

M = getManifold(Pose2)
np0 = getBelief(fg, :x0, :default) |> mean
np1 = getBelief(fg, :x1, :default) |> mean
np2 = getBelief(fg, :x2, :default) |> mean

# see issue JuliaRobotics/IncrementalInference.jl#1760
@test_broken isapprox(M, np0, ArrayPartition([0, 0.0], [0 -1.0; 1.0 0]), atol=1e-2)
@test_broken isapprox(M, np1, ArrayPartition([0, 1.0], [0 -1.0; 1.0 0]), atol=1e-2)
@test_broken isapprox(M, np2, ArrayPartition([0, 2.0], [0 -1.0; 1.0 0]), atol=1e-2)

##
end
