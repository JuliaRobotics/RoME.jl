using Distributed
using Dates
addprocs(4)
using RoME
using RoMEPlotting
using Gadfly
@everywhere using RoME

# Parse the arguments.
initial_timestep = parse(Int, ARGS[1])
final_timestep = parse(Int, ARGS[2])

# Let's load the Manhattan scenario using the g2o file.
file = (normpath(Base.find_package("RoME"), "../..", "examples", "manhattan_incremental.g2o"))
global instructions = importG2o(file)

# Make sure plots look a bit nicer.
latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt)
Gadfly.push_theme(latex_fonts)

function runManhattanStep(step::Integer)
    # Don't parse the instructions every time. Cache in global.
    global instructions

    # Choose where to save the step's data.
    data_logpath = "/tmp/caesar/tonio_results/$(step)-$(now())"
    # Create initial factor graph with specified logging path.
    fg = LocalDFG{SolverParams}(solverParams=SolverParams(logpath=data_logpath))

    # Add initial variable with a prior measurement to anchor the graph.
    addVariable!(fg, :x0, Pose2)
    initial_pose = MvNormal([0.0; 0.0; 0.0], Matrix(Diagonal([0.1;0.1;0.05].^2)))
    addFactor!(fg, [:x0], PriorPose2(initial_pose))

    # Add the desired amount of measurements given the step number.
    for i in 1:step
        parseG2oInstruction!(fg, instructions[i])
    end
    # And store a picture of the hitherto graph.
    padded_step = lpad(step, 4, "0")
    drawGraph(fg, show=false, engine="sfdp",
              filepath="$(getLogPath(fg))/graph$(padded_step).pdf")

    # Solve the graph, and save a copy of the tree.
    saveDFG(fg, "$(getLogPath(fg))/fg-before-solve")
    tree = solveTree!(fg)
    saveDFG(fg, "$(getLogPath(fg))/fg-after-solve")
    saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
    drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

    # Just store some quick plots.
    pl1 = drawPoses(fg, spscale=0.6)
    Gadfly.draw(PDF("$(getLogPath(fg))/poses$(padded_step).pdf", 20cm, 10cm), pl1)

    plkde = plotKDE(fg, ls(fg), dims=[1;2], levels=3)
    Gadfly.draw(PDF("$(getLogPath(fg))/kde$(padded_step).pdf", 20cm, 10cm), plkde)
end

# Run the loop using the command line arguments.
for i in initial_timestep:final_timestep
    runManhattanStep(i)
end
