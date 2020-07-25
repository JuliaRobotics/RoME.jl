using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using PackageCompiler

cd(@__DIR__)

create_sysimage(:RoME, sysimage_path="RoMESysimage.so", precompile_execution_file="precompile_triggers.jl")


## to use RoME with the newly created sysimage, start julia with:
# julia -J RoMESysimage.so
