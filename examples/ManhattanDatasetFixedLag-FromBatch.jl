using Distributed
using Dates
# addprocs(6) # Call instead with: julia -O3 -p18
using RoME
using RoMEPlotting
using Gadfly
@everywhere using RoME, RoMEPlotting, Gadfly

# Make sure plots look a bit nicer.
latex_fonts = Theme(major_label_font="CMU Serif", major_label_font_size=16pt,
                    minor_label_font="CMU Serif", minor_label_font_size=14pt,
                    key_title_font="CMU Serif", key_title_font_size=12pt,
                    key_label_font="CMU Serif", key_label_font_size=10pt)
Gadfly.push_theme(latex_fonts)

# Parse the arguments.
qfl_length = parse(Int, ARGS[1])

# Let's load the Manhattan scenario using the g2o file.
file = (normpath(Base.find_package("RoME"), "../..", "examples", "manhattan_incremental.g2o"))
global instructions = importG2o(file)

# Let's also get the file for batch factor graph at step 500.
fg_name = "fg-after-solve.tar.gz"
global fg_file = (normpath(Base.find_package("RoME"), "../..", "examples", fg_name))

function go_fixedlag_frombatch(qfl_length_arg::Integer)
    # Choose where to save the step's data.
    qfl_length = qfl_length_arg # Fixed lag window size.
    data_logpath = ENV["HOME"]*"/Documents/wafr/manhattan-frombatch-b$(qfl_length)-$(now())"

    # Instead of creating new graph, load the batch one.
    fg = LocalDFG{SolverParams}(solverParams=SolverParams(logpath=data_logpath))
    loadDFG(fg_file, Main, fg)
    tree = emptyBayesTree()

    # Set up the fixed lag smoothing.
    getSolverParams(fg).logpath = data_logpath
    getSolverParams(fg).isfixedlag = true
    getSolverParams(fg).qfl = qfl_length
    getSolverParams(fg).limitfixeddown = true
    getSolverParams(fg).dbg = true
    # CHECK I think we need to manually set all variables to not solvable.
    for vsym in ls(fg)
        getVariable(fg, vsym).solvable = 0
    end
    # But release the last 10 poses.
    activevars = []
    for i in 350:360
        push!(activevars, Symbol("x", i))
        getVariable(fg, Symbol("x", i)).solvable = 1
    end
    # CHECK As well as all of the factors.
    for fsym in lsf(fg)
        getFactor(fg, fsym).solvable = 0
    end
    # And activate the factors corresponding to the active variables.
    sym_act_vars = Array{Symbol, 1}(activevars)
    for fsym in getFactorsAmongVariablesOnly(fg, sym_act_vars)
        getFactor(fg, fsym).solvable = 1
    end

    # Add the next---or initial offset of---measurements to the graph.
    padded_step = lpad(501, 4, "0")
    parseG2oInstruction!(fg, instructions[501])

    # Solve the graph, and save a copy of the tree.
    saveDFG(fg, "$(getLogPath(fg))/fg-before-solve$(padded_step)")
    getSolverParams(fg).maxincidence = 1000
    tree = solveTree!(fg)
    saveDFG(fg, "$(getLogPath(fg))/fg-after-solve$(padded_step)")
    saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
    drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

    # Analyze clique counts.
    fid = open("$(getLogPath(fg))/clique-counts.txt", "w")
    nCliqs, nMarg, nReused, nBoth = calcCliquesRecycled(tree)
    println(fid, "$(padded_step), $(nCliqs), $(nMarg), $(nReused), $(nBoth)")

    # Just store some quick plots.
    pl1 = drawPoses(fg, spscale=0.6)
    Gadfly.draw(PDF("$(getLogPath(fg))/poses$(padded_step).pdf", 20cm, 10cm), pl1)

    # Solver stride.
    solveStride = 0
    # Run the loop for the remaining time steps.
    for i in 502:5453
        # Add the next measurement to the graph.
        parseG2oInstruction!(fg, instructions[i])
        padded_step = lpad(i, 4, "0")

        # Store each graph.
        saveDFG(fg, "$(getLogPath(fg))/fg-before-solve$(padded_step)")

        # Just store some quick plots, on another process
        remotecall((fgl, padded_stepl) -> begin
          @info "drawPoses, $(padded_stepl), for fg num variables=$(length(ls(fgl)))."
          pl1 = plotSLAM2DPoses(fgl, dyadScale=0.6, lbls=false)
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
        tree = solveTree!(fg, tree)
        saveDFG(fg, "$(getLogPath(fg))/fg-after-solve$(padded_step)")
        saveTree(tree, "$(getLogPath(fg))/tree$(padded_step).jld2")
        drawTree(tree, show=false, filepath="$(getLogPath(fg))/bt$(padded_step).pdf")

        # Analyze clique number.
        nCliqs, nMarg, nReused, nBoth = calcCliquesRecycled(tree)
        println(fid, "$(padded_step), $(nCliqs), $(nMarg), $(nReused), $(nBoth)")
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
go_fixedlag_frombatch(qfl_length)
