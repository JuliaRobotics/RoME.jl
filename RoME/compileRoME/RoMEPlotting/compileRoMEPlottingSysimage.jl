using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using PackageCompiler

cd(@__DIR__)

create_sysimage([:RoME,:RoMEPlotting], sysimage_path="RoMEPlottingSysimage.so", precompile_execution_file="precompile_triggers_plotting.jl")


## to use RoME and RoMEPlotting with the newly created sysimage, start julia with:
# julia -J RoMEPlottingSysimage.so
# February 2021: starting julia with the sysimage crashes
