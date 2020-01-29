using Distributed
addprocs(20)
using RoME
@everywhere using RoME

# Let's load the Manhattan scenario using the g2o file.
file = (normpath(Base.find_package("RoME"), "../..", "examples", "manhattan.g2o"))
instructions = importG2o(file)

# We can batch build the factor graph, or parse it line by line to recreate an
# incremental operation setting.
fg = LightDFG{SolverParams}(params=SolverParams())
num_of_measurements = 100
for i in 1:num_of_measurements
    parseG2oInstruction!(fg, instructions[i])
end
drawGraph(fg, show=true, engine="sfdp")

tree, smt, hist = solveTree!(fg)
