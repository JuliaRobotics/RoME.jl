# add more julia processes
# nprocs() < 3 ? addprocs(4-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using RoME, Distributions
## Inter-operating visualization packages for Caesar/RoME/IncrementalInference exist
using RoMEPlotting
using Compose
# Using this for plotting.
using Gadfly

# Standard Hexagonal example for totalIterations - solve every iterationsPerSolve iterations.
function runHexagonalExample(fg::FactorGraph, totalIterations::Int, iterationsPerSolve::Int)::Vector{Int}
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
            t = @timed batchSolve!(fg)
            push!(solveTimes, t[2])
        end
    end
    return solveTimes
end

# start with an empty factor graph object
fg = initfg()
# DO NOT enable fixed-lag operation
solveTimes = runHexagonalExample(fg, 100, 6)
# Plot the many iterations to see that it succeeded.
drawPosesLandms(fg)

# Plot the time taken per solve for the whole graph.
# This demonstrates batch solving linearity. We should only use this
# for offline, batch solving, and this forms the problem statement for
# why fixed-lag solving.
using Gadfly
Gadfly.plot(x=1:length(solveTimes), y=solveTimes, Geom.path,
    Guide.title("Solving Time vs. Iteration for Complete Graph"),
    Guide.xlabel("Solving Iteration"),
    Guide.ylabel("Solving Time (seconds)"))

#############################
##### ---- Fixed lag mode ---
#############################

fgFixedLag = initfg()
fgFixedLag.isfixedlag = true
fgFixedLag.qfl = 10

solveTimesFixedLag = runHexagonalExample(fgFixedLag, 100, 6)
# Plot the many iterations to see that it succeeded.
drawPosesLandms(fg)

# Plot the time taken per solve for the whole graph.
# This demonstrates batch solving linearity. We should only use this
# for offline, batch solving, and this forms the problem statement for
# why fixed-lag solving.
Gadfly.plot(x=1:length(solveTimesFixedLag), y=solveTimesFixedLag, Geom.path,
    Guide.title("Solving Time vs. Iteration for Fixed-Lag Operation"),
    Guide.xlabel("Solving Iteration"),
    Guide.ylabel("Solving Time (seconds)"))
