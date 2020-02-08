using RoME
using RoMEPlotting
using IncrementalInference
using DistributedFactorGraphs

# Find the filezzzzz.
fg_dir = "/home/tonio/repos/bayes-garden/data/batch-liljon-1500"
fg_name = "fg-after-solve.tar.gz"
fg_file = joinpath(fg_dir, fg_name)

# Load the factor graph.
fg = initfg()
loadDFG(fg_file, Main, fg)

# Draw some pretty plots.
pl1 = drawPoses(fg, spscale=0.6, lbls=false)
pl1 |> PDF(joinpath(fg_dir, "poses.pdf"), 20cm, 10cm)

# Useful functions to deal with sorting symbols.

# For crying out loud, find a way to lexicographically sort the symbols.
function natural(x, y)
    k(x) = [occursin(r"\d+", s) ? parse(Int, s) : s
            for s in split(replace(x, r"\d+" => s->" $s "))]
    A = k(x); B= k(y)
    for (a, b) in zip(A, B)
        if !isequal(a, b)
            return typeof(a) <: typeof(b) ? isless(a, b) :
                   isa(a,Int) ? true : false
        end
    end
    return length(A) < length(B)
end

function getSortedSymbols(fg::AbstractDFG, lexsortfxn)
    # First transform the factor graph symbols into strings to sort them
    # lexicographically, in case numbers are not padded.
    stringsyms = []
    for sym in sort(ls(fg))
        push!(stringsyms, "$(sym)")
    end
    sort!(stringsyms, lt=lexsortfxn)

    # Now transform sorted string array into sorted symbols array.
    sortedsyms = []
    for sym in stringsyms
        push!(sortedsyms, Symbol(sym))
    end
    return sortedsyms
end


# Output ground truth to file.
sortedsyms = getSortedSymbols(fg, natural)
fid = open(joinpath(fg_dir, "ground-truth.txt"), "w")
for i in 1:length(sortedsyms)
    @show i
    @show sortedsyms[i]
    @show getVariablePPE(fg, sortedsyms[i]).suggested
    x_hat = getVariablePPE(fg, sortedsyms[i]).suggested
    quat = normalize([0; 0; sin(x_hat[3]/2); cos(x_hat[3]/2)])
    println(fid, "$(i-1) $(x_hat[1]) $(x_hat[2]) 0.0 $(quat[1]) $(quat[2]) $(quat[3]) $(quat[4])")
end
close(fid)


# Encapsulate all of this into a function:
# Where `lexsort` is a function to sort lexicographically (e.g., `natural`).
function exportTrajEstimatesToFile(fg_path::String, lexsort)
    # Load the factor graph.
    fg = initfg()
    loadDFG(fg_path, Main, fg)

    # Draw some pretty plots.
    pl1 = drawPoses(fg, spscale=0.6, lbls=false)
    pl1 |> PDF(joinpath(dirname(fg_path), "trajectory-estimates-$(basename(dirname(fg_path))).pdf"), 20cm, 10cm)

    # Output ground truth to file.
    sortedsyms = getSortedSymbols(fg, lexsort)
    fid = open(joinpath(dirname(fg_path), "trajectory-estimates-$(basename(dirname(fg_path))).txt"), "w")
    for i in 1:length(sortedsyms)
        x_hat = getVariablePPE(fg, sortedsyms[i]).suggested
        quat = normalize([0; 0; sin(x_hat[3]/2); cos(x_hat[3]/2)])
        println(fid, "$(i-1) $(x_hat[1]) $(x_hat[2]) 0.0 $(quat[1]) $(quat[2]) $(quat[3]) $(quat[4])")
    end
    close(fid)
end

# Now run this for all of our runs:
b100_new_liljon_1411 = "/home/tonio/repos/bayes-garden/data/manhattan/b100-new-liljon-1411/fg-after-solve1411.tar.gz"
exportTrajEstimatesToFile(b100_new_liljon_1411, natural)

b10_old_augustus_1571 = "/home/tonio/repos/bayes-garden/data/manhattan/b10-old-augustus-1571/fg-after-solve1571.tar.gz"
exportTrajEstimatesToFile(b10_old_augustus_1571, natural)

b30_old_liljon_1451 = "/home/tonio/repos/bayes-garden/data/manhattan/b30-old-liljon-1451/fg-after-solve1451.tar.gz"
exportTrajEstimatesToFile(b30_old_liljon_1451, natural)

incremental_liljon_1101 = "/home/tonio/repos/bayes-garden/data/manhattan/incremental-liljon-1101/fg-after-solve1101.tar.gz"
exportTrajEstimatesToFile(incremental_liljon_1101, natural)

batch_liljon_1500 = "/home/tonio/repos/bayes-garden/data/manhattan/batch-liljon-1500/fg-after-solve.tar.gz"
exportTrajEstimatesToFile(batch_liljon_1500, natural)

batch_liljon_1800 = "/home/tonio/repos/bayes-garden/data/manhattan/batch-liljon-1800/fg-after-solve.tar.gz"
exportTrajEstimatesToFile(batch_liljon_1800, natural)
