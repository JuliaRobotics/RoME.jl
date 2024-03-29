using Distributed
using Dates
# addprocs(18) # call instead with:  julia -O3 -p18
using RoME
using RoMEPlotting
using Gadfly
@everywhere using RoME, RoMEPlotting, Gadfly

# Parse the arguments.
initial_offset = parse(Int, ARGS[1])
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

function go(initial_offset::Integer, final_timestep::Integer)
    # Choose where to save the step's data.
    data_logpath = "/media/data2/tonio_results/manhattan-$(now())"
    # Create initial factor graph with specified logging path.
    fg = LocalDFG{SolverParams}(solverParams=SolverParams(logpath=data_logpath))

    # adding debug
    getSolverParams(fg).dbg = true # store cliqSubFg at critical points during solve.

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

    # And store a picture of the hitherto graph.
    # drawGraph(fg, show=false, engine="sfdp",
    #           filepath="$(getLogPath(fg))/graph$(padded_step).pdf")

    # Solve the graph, and save a copy of the tree.
    saveDFG(fg, "$(getLogPath(fg))/fg-before-solve$(padded_step)")
    tree = solveTree!(fg)
    saveDFG(fg, "$(getLogPath(fg))/fg-after-solve$(padded_step)")
    saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
    drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

    # analyze cliques recycle
    fid = open("$(getLogPath(fg))/clique-counts.txt", "w")
    nMarg, nReused = calcCliquesRecycled(tree)
    nCliqs = length(tree.cliques)
    println(fid, "$(padded_step), $(nCliqs), $(nMarg), $(nReused)")
    # close fid at end

    # Just store some quick plots.
    pl1 = drawPoses(fg, spscale=0.6)
    Gadfly.draw(PDF("$(getLogPath(fg))/poses$(padded_step).pdf", 20cm, 10cm), pl1)

    # plkde = plotKDE(fg, ls(fg), dims=[1;2], levels=3)
    # Gadfly.draw(PDF("$(getLogPath(fg))/kde$(padded_step).pdf", 20cm, 10cm), plkde)

    # solve stride
    solveStride = 0
    # Run the loop for the remaining time steps.
    for i in (initial_offset + 1):final_timestep
        # Add the next measurement to the graph.
        parseG2oInstruction!(fg, instructions[i])


        padded_step = lpad(i, 4, "0")
        # And store a picture of the hitherto graph.
        # drawGraph(fg, show=false, engine="sfdp",
        #           filepath="$(getLogPath(fg))/graph$(padded_step).pdf")

        # store each graph
        saveDFG(fg, "$(getLogPath(fg))/fg-before-solve$(padded_step)")

        # Just store some quick plots, on another process
        remotecall((fgl, padded_stepl) -> begin
          @info "drawPoses, $(padded_stepl), for fg num variables=$(length(ls(fgl)))."
          pl1 = drawPoses(fgl, spscale=0.6, lbls=false)
          pl1 |> PDF("$(getLogPath(fgl))/poses$(padded_stepl).pdf", 20cm, 10cm)
        end, rand(Categorical(nprocs()-1))+1, fg, padded_step)

        # only solve every 10th instruction
        solveStride += 1
        if solveStride % 10 != 0
          @info "solveStride=$solveStride"
          continue
        end
        @info "Going for solve"

        # Solve the graph, and save a copy of the tree.
        getSolverParams(fg).maxincidence = 1000
        tree = solveTree!(fg, tree)
        saveDFG(fg, "$(getLogPath(fg))/fg-after-solve$(padded_step)")
        saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
        drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

        # Analyze clique number.
        nMarg, nReused = calcCliquesRecycled(tree)
        nCliqs = length(tree.cliques)
        println(fid, "$(padded_step), $(nCliqs), $(nMarg), $(nReused)")
        flush(fid)

        # plkde = plotKDE(fg, ls(fg), dims=[1;2], levels=3)
        # Gadfly.draw(PDF("$(getLogPath(fg))/kde$(padded_step).pdf", 20cm, 10cm), plkde)

        # force garbage collection to reduce memory footprint
        GC.gc()
    end
    # final plot
    padded_step = lpad(final_timestep+1, 4, "0")
    pl1 = drawPoses(fg, spscale=0.6)
    pl1 |> PDF("$(getLogPath(fg))/poses$(padded_step).pdf", 20cm, 10cm)

    close(fid)
end

# Run within a function to avoid undefined variable errors, and faster.
go(initial_offset, final_timestep)
