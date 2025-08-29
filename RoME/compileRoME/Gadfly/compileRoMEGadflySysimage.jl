using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using PackageCompiler

cd(@__DIR__)

create_sysimage([:RoME,:Gadfly,:GraphPlot,:DistributedFactorGraphs,:Compose], sysimage_path="RoMEGadflySysimage.so", precompile_execution_file="precompile_triggers_gadfly.jl")


## to use RoME Gadfly and the associated Pkgs with the newly created sysimage, start julia with:
# julia -J RoMEGadflySysimage.so
