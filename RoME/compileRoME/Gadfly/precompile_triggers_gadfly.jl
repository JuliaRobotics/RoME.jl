using GraphPlot, Compose #, ImageMagick, LightGraphs # RoME, DistributedFactorGraphs, Gadfly
# include(joinpath(pkgdir(RoME), "test", "runtests.jl")) # February 2021: some tests fail here
# include(joinpath(pkgdir(Gadfly), "test", "runtests.jl")) # February 2021: some tests fail here
include(joinpath(pkgdir(GraphPlot), "test", "runtests.jl"))
# include(joinpath(pkgdir(DistributedFactorGraphs), "test", "runtests.jl")) # February 2021: CI already errors?
# include(joinpath(pkgdir(Compose), "test", "runtests.jl")) # February 2021: some tests fail here
