
using Test
using Statistics
using RoME

##

@testset "Testing tree init prior usage" begin

##

fg = generateCanonicalFG_Circle(4;graphinit=true)
# deleteVariable!.(fg, [:x3, :x4, :l1])


# fg.solverParams.showtree = true
# fg.solverParams.drawtree = true
# fg.solverParams.dbg = true
fg.solverParams.useMsgLikelihoods = true

# now move the graph prior somewhere else (ie all init is wrong)
deleteFactor!(fg, :x0f1)
prpo = PriorPose2(MvNormal([5.,0.,0.0], 0.01*diagm([1;1;1.])))
addFactor!(fg, [:x0], prpo)

##

smtasks = Task[]
tree = solveTree!(fg, smtasks=smtasks); #, recordcliqs=ls(fg));

# hists = fetchCliqHistoryAll!(smtasks)
# printCSMHistoryLogical(hists)

##

X4 = getPoints(fg, :x4)
μX4 =  mean(getManifold(Pose2), X4)

@error("Must first fix IIF #913")
@test 2.0 < μX4.parts[1][1]
@test -1.5 < μX4.parts[1][2] < 1.5

##
end