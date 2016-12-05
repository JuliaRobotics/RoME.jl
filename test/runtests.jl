using RoME
using Base.Test


println("[TEST] Camera functions evaluations...")
include("testCameraFunctions.jl")
println("[SUCCESS]")

println("[TEST] Linear array functions evaluations...")
include("testDidsonFunctions.jl")
println("[SUCCESS]")
