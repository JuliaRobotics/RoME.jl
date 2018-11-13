# add more julia processes
using Distributed
nprocs() < 2 ? addprocs(2-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using RoME, Distributions
## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
# using RoMEPlotting
using Compose
# Using this for trend plotting.
using Gadfly
# Using for capturing data.
using DataFrames


numVariables = 300
solveEveryNVariables = 20
lagLength = 30

# Standard Hexagonal example for totalIterations - solve every iterationsPerSolve iterations.
function runHexagonalExample(fg::FactorGraph, totalIterations::Int, iterationsPerSolve::Int)::DataFrame
    # Add the first pose :x0
    addNode!(fg, :x0, Pose2)

    # Add at a fixed location PriorPose2 to pin :x0 to a starting location
    addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I, 3,3))))

    # Add a landmark l1
    addNode!(fg, :l1, Point2, labels=["LANDMARK"])

    # Drive around in a hexagon a number of times
    solveTimes = DataFrame(GraphSize = [], TimeBuildBayesTree = [], TimeSolveGraph = [])
    for i in 0:totalIterations
        psym = Symbol("x$i")
        nsym = Symbol("x$(i+1)")
        @info "Adding pose $nsym..."
        addNode!(fg, nsym, Pose2)
        pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal( [0.1;0.1;0.1].^2 ) )))
        @info "Adding odometry factor between $psym -> $nsym..."
        addFactor!(fg, [psym;nsym], pp )

        if i % 6 == 0
            @info "Creating factor between $psym and l1..."
            p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
            addFactor!(fg, [psym; :l1], p2br)
        end
        if i % iterationsPerSolve == 0 && i != 0
            @info "Performing inference!"
            if fg.isfixedlag
                @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
                fifoFreeze!(fg)
            end
            tBuild = @timed tree = wipeBuildNewTree!(fg)
            tInfer = @timed inferOverTree!(fg, tree, N=100)
            graphSize = length([ls(fg)[1]..., ls(fg)[2]...])
            push!(solveTimes, (graphSize, tBuild[2], tInfer[2]))
        end
    end
    return solveTimes
end

# start with an empty factor graph object
fg = initfg()
# DO NOT enable fixed-lag operation
solverTimesForBatch = runHexagonalExample(fg, numVariables, solveEveryNVariables)
# savejld(fg, "batchSolve.jld2")

#############################
##### ---- Fixed lag mode ---
#############################

fgFixedLag = initfg()
fgFixedLag.isfixedlag = true
fgFixedLag.qfl = lagLength

solverTimesFixedLag = runHexagonalExample(fgFixedLag, numVariables, solveEveryNVariables)
# savejld(fgFixedLag, "fixedSolve.jld2")

#### Visualization

# Plot the many iterations to see that it succeeded.
# Batch
# drawPosesLandms(fg)

# Fixed lag
# drawPosesLandms(fgFixedLag)

#### PLOTTING of Results

using Gadfly
using Colors
using CSV

# Make a clean dataset
rename!(solverTimesForBatch, :TimeBuildBayesTree => :Batch_BayedBuild, :TimeSolveGraph => :Batch_SolveGraph);
rename!(solverTimesFixedLag, :TimeBuildBayesTree => :FixedLag_BayedBuild, :TimeSolveGraph => :FixedLag_SolveGraph);
timingMerged = DataFrames.join(solverTimesForBatch, solverTimesFixedLag, on=:GraphSize)
CSV.write("timing_comparison.csv", timingMerged)

# PP = []
# push!(PP, Gadfly.layer(x=timingMerged[:GraphSize], y=timingMerged[:FixedLag_SolveGraph], Geom.path, Theme(default_color=colorant"green"))[1]);
# push!(PP, Gadfly.layer(x=timingMerged[:GraphSize], y=timingMerged[:Batch_SolveGraph], Geom.path, Theme(default_color=colorant"magenta"))[1]);
#
# plt = Gadfly.plot(PP...,
#     Guide.title("Solving Time vs. Iteration for Fixed-Lag Operation"),
#     Guide.xlabel("Solving Iteration"),
#     Guide.ylabel("Solving Time (seconds)"),
#     Guide.manual_color_key("Legend", ["fixed", "batch"], ["green", "magenta"]))
# Gadfly.draw(PNG("results_comparison.png", 12cm, 15cm), plt)


#### EXPORT the data
using JLD2, FileIO
@save "results.jld2" timingMerged
