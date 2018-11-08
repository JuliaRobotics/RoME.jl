# add more julia processes
# nprocs() < 3 ? addprocs(4-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using RoME, Distributions
## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
# using RoMEPlotting
using Compose
# Using this for trend plotting.
using Gadfly


numVariables = 200
solveEveryNVariables = 6
lagLength = 15

# Standard Hexagonal example for totalIterations - solve every iterationsPerSolve iterations.
function runHexagonalExample(fg::FactorGraph, totalIterations::Int, iterationsPerSolve::Int)#::Vector{Float64}
    # Add the first pose :x0
    addNode!(fg, :x0, Pose2)

    # Add at a fixed location PriorPose2 to pin :x0 to a starting location
    addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I, 3,3))))

    # Add a landmark l1
    addNode!(fg, :l1, Point2, labels=["LANDMARK"])

    # Drive around in a hexagon a number of times
    solveTimes = []
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
        if i % iterationsPerSolve == 0
            @info "Performing inference!"
            if fg.isfixedlag
                @info "Quasi fixed-lag is enabled (a feature currently in testing)!"
                fifoFreeze!(fg)
            end
            tBuild = @timed tree = wipeBuildNewTree!(fg)
            tInfer = @timed inferOverTree!(fg, tree, N=100)
            push!(solveTimes, [tBuild[2], tInfer[2]])
        end
    end
    return solveTimes
end

# start with an empty factor graph object
fg = initfg()
# DO NOT enable fixed-lag operation
solverTimesForBatch = runHexagonalExample(fg, numVariables, solveEveryNVariables)
# t = @timed batchSolve!(fg)
# Plot the many iterations to see that it succeeded.
# drawPosesLandms(fg)

# Plot the time taken per solve for the whole graph.
# This demonstrates batch solving linearity. We should only use this
# for offline, batch solving, and this forms the problem statement for
# why fixed-lag solving.
# using Gadfly
# Gadfly.plot(x=1:length(solveTimes), y=solveTimes, Geom.path,
#     Guide.title("Solving Time vs. Iteration for Complete Graph"),
#     Guide.xlabel("Solving Iteration"),
#     Guide.ylabel("Solving Time (seconds)"))

#############################
##### ---- Fixed lag mode ---
#############################

# Batch control for 60 - no fixed lag and [buildBayes, infer]
# solverTimesForBatch = [[0.339419, 2.92293],
# [1.20682, 9.57546],
# [2.13343, 17.3296],
# [3.06916, 26.7324],
# [4.0266, 37.0113],
# [5.12491, 47.7176],
# [7.0712, 75.0297],
# [10.116, 88.8075],
# [11.8314, 104.447],
# [13.2786, 118.219],
# [14.3143, 127.049]]

# Fixed lag
# [5.0763, 13.0591]
# [0.0895069, 6.71251]
# [0.0607179, 8.96031]
# [0.140903, 10.9674]
# [0.311809, 12.6977]
# [0.669601, 15.1875]
# [1.17332, 17.8367]
# [1.78019, 20.373]
# [2.63709, 24.3739]
# [3.773, 29.1994]
# [5.10071, 34.9171]

fgFixedLag = initfg()
fgFixedLag.isfixedlag = true
fgFixedLag.qfl = lagLength

solveTimesFixedLag = runHexagonalExample(fgFixedLag, numVariables, solveEveryNVariables)

# tFixed = @timed batchSolve!(fgFixedLag)
# @show tFixed[2]
# Plot the many iterations to see that it succeeded.
# drawPosesLandms(fg)

# Plot the time taken per solve for the whole graph.
# This demonstrates batch solving linearity. We should only use this
# for offline, batch solving, and this forms the problem statement for
# why fixed-lag solving.
# Gadfly.plot(x=1:length(inferTimes), y=inferTimes, Geom.path,
#     Guide.title("Solving Time vs. Iteration for Fixed-Lag Operation"),
#     Guide.xlabel("Solving Iteration"),
#     Guide.ylabel("Solving Time (seconds)"))

#### PLOTTING

using Gadfly

using Colors

PP = []

batchTimes = [solverTimesForBatch[i][2] for i in 1:length(solverTimesForBatch)]
inferTimes = [solveTimesFixedLag[i][2] for i in 1:length(solveTimesFixedLag)]

push!(PP, Gadfly.layer(x=1:length(inferTimes), y=inferTimes, Geom.path, Theme(default_color=colorant"green"))[1]);
push!(PP, Gadfly.layer(x=1:length(batchTimes), y=batchTimes, Geom.path, Theme(default_color=colorant"magenta"))[1]);

Gadfly.plot(PP...,
    Guide.title("Solving Time vs. Iteration for Fixed-Lag Operation"),
    Guide.xlabel("Solving Iteration"),
    Guide.ylabel("Solving Time (seconds)"),
    Guide.manual_color_key("Legend", ["fixed", "batch"], ["green", "magenta"]))
