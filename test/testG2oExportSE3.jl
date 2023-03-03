
# using Revise

using RoME
using OrderedCollections: OrderedDict
using Test

##

@testset "Test g2o SE3 export" begin
##

# fg = initfg()
# addVariable!.(fg, [:x0;:x1;:x2;:x3], Pose3)
# addFactor!(fg, [:x0], PriorPose3(MvNormal(zeros(6),diagm(0.1*ones(6)))))
# addFactor!(fg, [:x0;:x1], Pose3Pose3(MvNormal([1;0;0;0;0;0.],diagm(0.1*ones(6)))); graphinit=false)
# addFactor!(fg, [:x1;:x2], Pose3Pose3(MvNormal([1;0;0;0;0;0.],diagm(0.1*ones(6)))); graphinit=false)
# addFactor!(fg, [:x2;:x3], Pose3Pose3(MvNormal([1;0;0;0;0;0.],diagm(0.1*ones(6)))); graphinit=false)
# saveDFG(joinpath(@__DIR__,"testdata","g2otest.tar.gz"), fg)

fg = loadDFG(joinpath(@__DIR__,"testdata","g2otest.tar.gz"))

setPPE!.(fg, ls(fg), :parametric)
g2ofile = joinpath("/tmp", "caesar", "export.g2o")
mkpath(dirname(g2ofile))
varIntLabel = OrderedDict(zip(ls(fg), (0:length(ls(fg))-1)))

##

g2ofile = "/tmp/test.g2o"
RoME.exportG2o(fg; filename=g2ofile, varIntLabel, solveKey = :parametric)

Base.rm(g2ofile)

##
end

##
