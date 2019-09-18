using Test
using DistributedFactorGraphs
using IncrementalInference, RoME

# Make a simple graph
dfg = GraphsDFG{SolverParams}(params=SolverParams())
# Add the first pose :x0
x0 = addVariable!(dfg, :x0, Pose2)
# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
prior = addFactor!(dfg, [:x0], PriorPose2( MvNormal([10; 10; 1.0/8.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) ) )
# Drive around in a hexagon
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(dfg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(dfg, [psym;nsym], pp )
end

# Save it
saveFolder = "/tmp/fileDFG"
saveDFG(dfg, saveFolder)
@test readdir("$saveFolder/variables") == ["x0.json", "x1.json", "x2.json", "x3.json", "x4.json", "x5.json", "x6.json"]
@test readdir("$saveFolder/factors") == ["x0f1.json", "x0x1f1.json", "x1x2f1.json", "x2x3f1.json", "x3x4f1.json", "x4x5f1.json", "x5x6f1.json"]

retDFG = GraphsDFG{SolverParams}(params=SolverParams())
retDFG = loadDFG(saveFolder, IncrementalInference, retDFG)
@test symdiff(ls(dfg), ls(dfg)) == []
@test symdiff(lsf(dfg), lsf(retDFG)) == []

# Solve it
tree, smt, hist = solveTree!(retDFG)

# Checking estimates exist for all variables.
for variable in getVariables(retDFG)
    @info "Testing if $(variable.label) has estimates... = $(haskey(variable.estimateDict, :default))"
    @test haskey(variable.estimateDict, :default)
    if haskey(variable.estimateDict, :default)
        @test haskey(variable.estimateDict[:default], :max)
        @test haskey(variable.estimateDict[:default], :mean)
        @test haskey(variable.estimateDict[:default], :ppe)
    end
end

# Save it and load it again
saveDFG(retDFG, saveFolder)
retDFG = GraphsDFG{SolverParams}(params=SolverParams())
retDFG = loadDFG(saveFolder, IncrementalInference, retDFG)

for variable in getVariables(retDFG)
    @info "Testing if $(variable.label) has estimates... = $(haskey(variable.estimateDict, :default))"
    @test haskey(variable.estimateDict, :default)
    if haskey(variable.estimateDict, :default)
        @test haskey(variable.estimateDict[:default], :max)
        @test haskey(variable.estimateDict[:default], :mean)
        @test haskey(variable.estimateDict[:default], :ppe)
    end
end

# Success!
