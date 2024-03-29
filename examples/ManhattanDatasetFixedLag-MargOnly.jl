using Distributed
using Dates
# addprocs(6) # Call instead with: julia -O3 -p18
using RoME
using RoMEPlotting
using Gadfly
@everywhere using RoME, RoMEPlotting, Gadfly

# Parse the arguments.
initial_offset = parse(Int, ARGS[1])
final_timestep = parse(Int, ARGS[2])
qfl_length = parse(Int, ARGS[3])

# Let's load the Manhattan scenario using the g2o file.
file = (normpath(Base.find_package("RoME"), "../..", "examples", "manhattan_incremental.g2o"))
global instructions = importG2o(file)

# Make sure plots look a bit nicer.
latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt)
Gadfly.push_theme(latex_fonts)

function go_fixedlag(initial_offset::Integer, final_timestep::Integer, qfl_length_arg::Integer)
    # Choose where to save the step's data.
    qfl_length = qfl_length_arg # Fixed lag window size.
    data_logpath = "/media/data2/tonio_results/manhattan-margONLY-qfl$(qfl_length)-$(now())"

    # Create initial factor graph with specified logging path.
    fg = LocalDFG{SolverParams}(solverParams=SolverParams(logpath=data_logpath))
    tree = emptyBayesTree()

    # Set up the fixed lag smoothing.
    getSolverParams(fg).isfixedlag = true
    getSolverParams(fg).qfl = qfl_length
    getSolverParams(fg).limitfixeddown = true
    getSolverParams(fg).dbg = true

    # Add initial variable with a prior measurement to anchor the graph.
    addVariable!(fg, :x0, Pose2)
    initial_pose = MvNormal([0.0; 0.0; 0.0], Matrix(Diagonal([0.1;0.1;0.05].^2)))
    addFactor!(fg, [:x0], PriorPose2(initial_pose))

    # Add the next---or initial offset of---measurements to the graph.
    padded_step = lpad(1, 4, "0")
    if initial_offset == 1
        parseG2oInstruction!(fg, instructions[1])
    else
        for j in 1:initial_offset
            parseG2oInstruction!(fg, instructions[j])
        end
        padded_step = lpad(initial_offset, 4, "0")
    end

    # Solve the graph, and save a copy of the tree.
    saveDFG(fg, "$(getLogPath(fg))/fg-before-solve$(padded_step)")
    getSolverParams(fg).maxincidence = 1000
    tree = solveTree!(fg)
    saveDFG(fg, "$(getLogPath(fg))/fg-after-solve$(padded_step)")
    saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
    drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

    # Analyze clique counts.
    fid = open("$(getLogPath(fg))/clique-counts.txt", "w")
    nMarg, nReused = calcCliquesRecycled(tree)
    nCliqs = length(tree.cliques)
    println(fid, "$(padded_step), $(nCliqs), $(nMarg), $(nReused)")

    # Just store some quick plots.
    pl1 = drawPoses(fg, spscale=0.6)
    Gadfly.draw(PDF("$(getLogPath(fg))/poses$(padded_step).pdf", 20cm, 10cm), pl1)

    # Solver stride.
    solveStride = 0
    # Run the loop for the remaining time steps.
    for i in (initial_offset + 1):final_timestep
        # Add the next measurement to the graph.
        parseG2oInstruction!(fg, instructions[i])
        padded_step = lpad(i, 4, "0")

        # Store each graph.
        saveDFG(fg, "$(getLogPath(fg))/fg-before-solve$(padded_step)")

        # Just store some quick plots, on another process
        remotecall((fgl, padded_stepl) -> begin
          @info "drawPoses, $(padded_stepl), for fg num variables=$(length(ls(fgl)))."
          pl1 = plotSLAM2DPoses(fgl, spscale=0.6, lbls=false)
          pl1 |> PDF("$(getLogPath(fgl))/poses$(padded_stepl).pdf", 20cm, 10cm)
        end, rand(Categorical(nprocs()-1))+1, fg, padded_step)

        # Only solve every 10th instruction.
        solveStride += 1
        if solveStride % 10 != 0
            @info "solveStride=$solveStride"
            continue
        end
        @info "Going for solve"

        # Solve the graph, and save a copy of the tree.
        getSolverParams(fg).maxincidence = 1000
        tree = solveTree!(fg)
        saveDFG(fg, "$(getLogPath(fg))/fg-after-solve$(padded_step)")
        saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
        drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

        # Analyze clique number.
        nMarg, nReused = calcCliquesRecycled(tree)
        nCliqs = length(tree.cliques)
        println(fid, "$(padded_step), $(nCliqs), $(nMarg), $(nReused)")
        flush(fid)

        # Run the garbage collector.
        GC.gc()
    end
    # Final plot.
    padded_step = lpad(final_timestep+1, 4, "0")
    pl1 = drawPoses(fg, spscale=0.6)
    pl1 |> PDF("$(getLogPath(fg))/poses$(padded_step).pdf", 20cm, 10cm)

    close(fid)
end

# Run within a function to avoid undefined variable errors.
go_fixedlag(initial_offset, final_timestep, qfl_length)
