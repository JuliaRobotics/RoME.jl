
# using Revise

using RoME
using OrderedCollections: OrderedDict
using Test

##

@testset "Test g2o SE3 export" begin
##

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
