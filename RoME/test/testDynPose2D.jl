using RoME
using Test
using TensorCast
using Dates
using LinearAlgebra
using DistributedFactorGraphs.JSON
using DistributedFactorGraphs

##

@testset "test DynPose2 and velocity..." begin

##

    N = 100
    fg = initfg()

    # add first pose locations
    addVariable!(fg, :x0, DynPose2; timestamp = DFG.Timestamp(Nanosecond(0)))

    # Prior factor as boundary condition
    pp0 = DynPose2VelocityPrior(
        MvNormal(zeros(3), Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([10.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    addFactor!(fg, [:x0;], pp0)

    # initialize the first pose
    IncrementalInference.doautoinit!(fg, [getVariable(fg, :x0);])

    addVariable!(fg, :x1, DynPose2; timestamp = DFG.Timestamp(Nanosecond(1000_000_000)))

    # conditional likelihood between Dynamic Point2
    dp2dp2 = VelPose2VelPose2(
        MvNormal([10.0; 0; 0], Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    addFactor!(fg, [:x0; :x1], dp2dp2)

    pts = approxConv(fg, :x0x1f1, :x1)

    # initAll!(fg)

    tree = solveTree!(fg)

    _X1 = DFG.getCoordinates.(DynPose2, getVal(fg, :x1))
    @cast X1[j, i] := _X1[i][j]
    @test 0.9 * N <= sum(abs.(X1[1, :] .- 10.0) .< 0.75)
    @test 0.9 * N <= sum(abs.(X1[2, :] .- 0.0) .< 0.75)
    # @show TU.wrapRad.(X1[3,:])
    @test 0.8 * N <= sum(abs.(TU.wrapRad.(X1[3, :]) .- 0.0) .< 0.25)
    # @warn "wrapRad issue, accepting 80% as good enough until issue JuliaRobotics/RoME.jl#90 is fixed."
    @test 0.9 * N <= sum(abs.(X1[4, :] .- 10.0) .< 0.5)
    @test 0.9 * N <= sum(abs.(X1[5, :] .- 0.0) .< 0.5)

    # using RoMEPlotting
    # # plotLocalProduct(fg, :x10, dims=[1;2])
    # # plotSLAM2DPoses(fg)
    # xx1 = marginal(getBelief(fg, :x1),[1;2;3])
    # plotPose(Pose2(), [xx1])
    # plotPose(fg, [:x1], levels=1, show=false)
    #
    # plotKDE(marginal(getBelief(fg, :x1),[4;5]), levels=5)

##

end

@testset "test distribution compare functions..." begin

##

    mu = randn(6)
    mv1 = MvNormal(deepcopy(mu), Matrix{Float64}(LinearAlgebra.I, 6, 6))
    mv2 = MvNormal(deepcopy(mu), Matrix{Float64}(LinearAlgebra.I, 6, 6))
    mv3 = MvNormal(randn(6), Matrix{Float64}(LinearAlgebra.I, 6, 6))
    @test RoME.compareDensity(mv1, mv2)
    @test !RoME.compareDensity(mv1, mv3)
    @test !RoME.compareDensity(mv2, mv3)

##

end

@testset "test DynPose2 packing converters..." begin

##

    pp0 = DynPose2VelocityPrior(
        MvNormal(zeros(3), Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([10.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    from_fac = FactorDFG((:x1,), pp0)
    jstr = JSON.json(from_fac; style = DFG.DFGJSONStyle())
    to_fac = JSON.parse(jstr, FactorDFG; style = DFG.DFGJSONStyle())
    # check if the original and unpacked are equal
    @test from_fac == to_fac

    #
    dp2dp2 = VelPose2VelPose2(
        MvNormal([10.0; 0; 0], Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    from_fac = FactorDFG((:x1, :x2), dp2dp2)
    jstr = JSON.json(from_fac; style = DFG.DFGJSONStyle())
    to_fac = JSON.parse(jstr, FactorDFG; style = DFG.DFGJSONStyle())
    # check if the original and unpacked are equal
    @test from_fac == to_fac

##

end

@testset "test many DynPose2 chain stationary and 'pulled'..." begin

##

    N = 100
    fg = initfg()

    # add first pose locations
    addVariable!(fg, :x0, DynPose2; timestamp = DFG.Timestamp(Nanosecond(0)))

    # Prior factor as boundary condition
    pp0 = DynPose2VelocityPrior(
        MvNormal(zeros(3), Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    addFactor!(fg, [:x0;], pp0)

    sym = :x0
    k = 0
    for sy in Symbol[Symbol("x$i") for i = 1:10]
        k += 1
        addVariable!(fg, sy, DynPose2; timestamp = DFG.Timestamp(Nanosecond(1000_000_000 * k)))

        # conditional likelihood between Dynamic Point2
        dp2dp2 = VelPose2VelPose2(
            MvNormal([0.0; 0; 0], Matrix(Diagonal([1.0; 0.1; 0.001] .^ 2))),
            MvNormal([0.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
        )
        addFactor!(fg, [sym; sy], dp2dp2)
        sym = sy
    end # for

    initAll!(fg)
    x5 = mean(getBelief(getVariable(fg, :x5)))

    @test abs(x5[1]) < 1.25
    @test abs(x5[2]) < 1.25
    @test_broken abs(TU.wrapRad(x5[3])) < 0.4
    @test abs(x5[4]) < 0.5
    @test abs(x5[5]) < 0.5

    initAll!(fg)

    x10 = mean(getBelief(getVariable(fg, :x10)))

    @test abs(x10[1]) < 1.25
    @test abs(x10[2]) < 1.25
    @test_broken abs(TU.wrapRad(x10[3])) < 0.4
    @test abs(x10[4]) < 0.5
    @test abs(x10[5]) < 0.5

    # drawGraph(fg, show=true)
    # tree = buildTreeReset!(fg)
    # drawTree(tree, show=true)

    # using RoMEPlotting
    # Gadfly.set_default_plot_size(35cm, 25cm)
    # plotSLAM2DPoses(fg)
    # plotPose(fg, [:x10])

##

    @error ".useMsgLikelihoods = false required until IIF #1010 completed."
    getSolverParams(fg).useMsgLikelihoods = false

    # solve
    smtasks = Task[]
    tree = solveTree!(fg; smtasks = smtasks) #, recordcliqs=ls(fg));

##

    x5 = IIF.calcMeanMaxSuggested(getVariable(fg, :x5)).suggested
    # x5 = KDE.getKDEMean(getBelief(getVariable(fg, :x5)))

    @test abs(x5[1]) < 1.5
    @test abs(x5[2]) < 1.5
    @test_broken abs(TU.wrapRad(x5[3])) < 0.4
    @test abs(x5[4]) < 0.5
    @test abs(x5[5]) < 0.5

    x10 = mean(getBelief(getVariable(fg, :x10)))

    @test abs(x10[1]) < 2.75
    @test abs(x10[2]) < 2.75
    @test_broken abs(TU.wrapRad(x10[3])) < 0.5
    @test abs(x10[4]) < 0.5
    @test abs(x10[5]) < 0.5

    # pull the tail end out with position
    pp10 = DynPose2VelocityPrior(
        MvNormal([10.0; 0; 0], Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    addFactor!(fg, [:x10;], pp10)

    fg2 = deepcopy(fg)

    tree = solveTree!(fg) # N=N

##

    x10 = mean(getBelief(getVariable(fg, :x10)))

    @test 5.0 < x10[1]
    @test abs(x10[2]) < 1.0
    @test_broken abs(TU.wrapRad(x10[3])) < 0.6
    @test -0.1 < x10[4] < 1.0
    @test abs(x10[5]) < 0.5

    for sym in [Symbol("x$i") for i = 2:9]
        XX = mean(getBelief(getVariable(fg, sym)))

        @show sym, round.(XX, digits = 5)
        @test -2.0 < XX[1] < 10.0
        @test abs(XX[2]) < 1.0
        @test abs(TU.wrapRad(XX[3])) < 1.3
        @test -0.5 < XX[4] < 2.0
        @test abs(XX[5]) < 0.5
    end

##

end

# using RoMEPlotting
# plotLocalProduct(fg, :x10, dims=[1;2])
# drawPoses(fg)
# plotPose(fg, [:x9],levels=1);
#
# plotPose(fg, [:x1;:x2;:x3;:x4;:x5;:x6;:x7;:x8;:x9;:x10],levels=1);

# savejld(fg) # tempfg.jld

@testset "test many DynPose2 sideways velocity..." begin
##

    N = 100
    fg = initfg()

    # add first pose locations
    addVariable!(fg, :x0, DynPose2; timestamp = DFG.Timestamp(Nanosecond(0)))

    # Prior factor as boundary condition
    pp0 = DynPose2VelocityPrior(
        MvNormal([0.0; 0.0; pi / 2], Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.5; 0.5] .^ 2))),
    )
    addFactor!(fg, [:x0;], pp0)

    addVariable!(fg, :x1, DynPose2; timestamp = DFG.Timestamp(Nanosecond(1000_000_000)))

    pp0 = DynPose2VelocityPrior(
        MvNormal([1.0; 0.0; pi / 2], Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.5; 0.5] .^ 2))),
    )
    addFactor!(fg, [:x1;], pp0)

    # conditional likelihood between Dynamic Point2
    dp2dp2 = VelPose2VelPose2(
        MvNormal([0.0; -1.0; 0], Matrix(Diagonal([0.01; 0.01; 0.001] .^ 2))),
        MvNormal([0.0; 0], Matrix(Diagonal([0.1; 0.1] .^ 2))),
    )
    addFactor!(fg, [:x0; :x1], dp2dp2)

    getSolverParams(fg).N = N
    solveTree!(fg)
##
    # test for velocity in the body frame
    x0 = mean(getBelief(getVariable(fg, :x0)))

    @test -0.4 < x0[1] < 2.0
    @test abs(x0[2]) < 0.5
    @test_broken abs(x0[3] - pi / 2) < 0.1
    @test_broken abs(x0[4]) < 0.4
    @test -1.5 < x0[5] < -0.5

    x1 = mean(getBelief(getVariable(fg, :x1)))

    @test -0.1 < x1[1] < 2.0
    @test abs(x1[2]) < 0.5
    @test_broken abs(x1[3] - pi / 2) < 0.1
    @test_broken abs(x1[4]) < 0.4
    @test -1.5 < x1[5] < -0.5

##

end

#
