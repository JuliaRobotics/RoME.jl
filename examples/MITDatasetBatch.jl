using Distributed
using Dates
using RoME
using RoMEPlotting
using Gadfly
@everywhere using RoME, RoMEPlotting, Gadfly

total_meas = parse(Int, ARGS[1])
# Max number is 827.

# Let's load the Manhattan scenario using the g2o file.
file = (normpath(Base.find_package("RoME"), "../..", "examples", "MIT_incremental.g2o"))
global instructions = importG2o(file)

# Make sure plots look a bit nicer.
latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt)
Gadfly.push_theme(latex_fonts)

function solve_batch(total_meas::Integer)
    # Choose where to save the data.
    mit_total_meas = total_meas
    data_logpath = ENV["HOME"]*"/Documents/wafr/mit-batch-$(total_meas)-$(now())"

    # Create initial factor graph with specified logging path.
    fg = LightDFG{SolverParams}(params=SolverParams(logpath=data_logpath))

    # Add initial variable with a prior measurement to anchor the graph.
    addVariable!(fg, :x0, Pose2)
    initial_pose = MvNormal([0.0; 0.0; 0.0], Matrix(Diagonal([0.1;0.1;0.05].^2)))
    addFactor!(fg, [:x0], PriorPose2(initial_pose))

    # Add all variables and measurements.
    for i in 1:mit_total_meas
        parseG2oInstruction!(fg, instructions[i])
    end

    # Solve the graph, and save a copy of the tree.
    saveDFG(fg, "$(getLogPath(fg))/fg-before-solve")
    tree, smt, hist = solveTree!(fg, maxparallel=1000)
    saveDFG(fg, "$(getLogPath(fg))/fg-after-solve")
    saveTree(tree, "$(getLogPath(fg))/tree$(mit_total_meas).jld2")
    drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt.pdf")

    # Just store some quick plots.
    pl1 = drawPoses(fg, spscale=0.6, lbls=false)
    Gadfly.draw(PDF("$(getLogPath(fg))/poses$(mit_total_meas).pdf", 20cm, 10cm), pl1)

    # Run the garbage collector.
    GC.gc()
end

solve_batch(total_meas)
