
using Dates
using RoME
using CSV

df = CSV.read(@__DIR__() * "/testdata/IMUCase1.csv", DataFrame)
σ_a = 5e-4 # noise density m/s²/√Hz
σ_ω = 2e-6  # noise density rad/√Hz
Σy  = diagm([ones(3)*σ_a^2; ones(3)*σ_ω^2])


a_b = SA[0.,0,0]
ω_b = SA[0.,0,0]

accels = map(SVector, zip(df.nax, df.nay, df.naz))
gyros = map(SVector, zip(df.nwx, df.nwy, df.nwz))
timestamps = df.time

pidx = range(1; step=100, stop=length(timestamps))
factors = map(zip(pidx[1:end-1], pidx[2:end])) do (fr,to)
    r = range(fr,to-1)
    RoME.IMUDeltaFactor(
        SVector{3,Float64}.(accels[r]),
        SVector{3,Float64}.(gyros[r]),
        ones(length(r))*Ts,
        Σy,
        a_b,
        ω_b
    )
end

times = timestamps[pidx]

fg = initfg()
fg.solverParams.graphinit = false

foreach(enumerate(Nanosecond.(times * 10^9))) do (i,nanosecondtime)
    addVariable!(fg, Symbol("x",i-1), Pose3; nanosecondtime)
    addVariable!(fg, Symbol("v",i-1), Position3; nanosecondtime)
end

addFactor!(fg, [:x0], PriorPose3(MvNormal([1.0, 1, 1, 0, 0, 0], diagm(ones(6)*1e-3^2))))
# addFactor!(fg, [:v0], PriorPoint3(MvNormal([0.0, 0, 0], diagm(ones(3)*1e-3))))

# addFactor!(fg, [:x9], PriorPose3(MvNormal([1.0, 1, 1.0, 0, 0, 0], diagm(ones(6)*1e-3^2))))
# addFactor!(fg, [:v9], PriorPoint3(MvNormal([0.0, 0, 0], diagm(ones(3)*1e-3))))

for l in ls(fg, r"v\d")
    addFactor!(fg, [l], PriorPoint3(MvNormal([0.0, 0, 0], diagm(ones(3)*1e-3))))
end

# Variable for bias 
addVariable!(fg, :b, Position{6})

for (i,fac) in enumerate(factors)
    frx = Symbol("x",i-1)
    frv = Symbol("v",i-1)
    tox = Symbol("x",i)
    tov = Symbol("v",i)
    if exists(fg, :b)
        addFactor!(fg, [frx,frv, tox, tov, :b], fac)
    else
        addFactor!(fg, [frx,frv, tox, tov], fac)
    end
end


debug = [:Iteration, :Change,  " | ", :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :damping_term, "\n", :Stop]
debug = [:Iteration, :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :last_step_successful, "\n", :Stop]
stopping_criterion=StopAfterIteration(1000) | StopWhenGradientNormLess(1e-12) | StopWhenStepsizeLess(1e-12)

@time IIF.solveGraphParametric!(fg; stopping_criterion, debug, is_sparse=false, damping_term_min=1e-12, expect_zero_residual=false);

vnd = getVariableSolverData(fg, :b, :parametric)
vnd.val

## Test Case 2
df = CSV.read(@__DIR__() * "/testdata/IMUCase2.csv", DataFrame; limit = 3001)
# df = CSV.read(@__DIR__() * "/testdata/IMUCase2.csv", DataFrame)
σ_a = 0.037/60*sqrt(100) # noise density m/s²/√Hz - 0.037 m/s/√hr
σ_ω = 0.15*pi*60/180/3600*sqrt(100)  # noise density rad/√Hz - 0.15 deg/√hr
Σy  = diagm([ones(3)*σ_a^2; ones(3)*σ_ω^2])

a_b = SA[0.,0,0]
ω_b = SA[0.,0,0]

accels = map(SVector, zip(df.nax, df.nay, df.naz))
gyros = map(SVector, zip(df.nwx, df.nwy, df.nwz))
timestamps = df.time

pidx = range(1; step=100, stop=length(timestamps))
factors = map(zip(pidx[1:end-1], pidx[2:end])) do (fr,to)
    r = range(fr,to-1)
    RoME.IMUDeltaFactor(
        SVector{3,Float64}.(accels[r]),
        SVector{3,Float64}.(gyros[r]),
        ones(length(r))*Ts,
        Σy,
        a_b,
        ω_b
    )
end

times = timestamps[pidx]

fg = initfg()
fg.solverParams.graphinit = false

foreach(enumerate(Nanosecond.(times * 10^9))) do (i,nanosecondtime)
    addVariable!(fg, Symbol("x",i-1), Pose3; nanosecondtime)
    addVariable!(fg, Symbol("v",i-1), Position3; nanosecondtime)
end

addFactor!(fg, [:x0], PriorPose3(MvNormal([1.0, 1, 1, 0, 0, 0], diagm(ones(6)*1e-3^2))))
addFactor!(fg, [sortDFG(ls(fg, r"^x"))[end]], PriorPose3(MvNormal([1.0, 1, 1.0, 0, 0, 0], diagm(ones(6)*1e-3^2))))

for l in ls(fg, r"v\d")
    addFactor!(fg, [l], PriorPoint3(MvNormal([0.0, 0, 0], diagm(ones(3)*1e-3))))
end

# Variable for bias 
addVariable!(fg, :b, Position{6})

for (i,fac) in enumerate(factors)
    frx = Symbol("x",i-1)
    frv = Symbol("v",i-1)
    tox = Symbol("x",i)
    tov = Symbol("v",i)
    if exists(fg, :b)
        addFactor!(fg, [frx,frv, tox, tov, :b], fac)
    else
        addFactor!(fg, [frx,frv, tox, tov], fac)
    end
end


debug = [:Iteration, :Change,  " | ", :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :damping_term, "\n", :Stop]
debug = [:Iteration, :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :last_step_successful, "\n", :Stop]
stopping_criterion=StopAfterIteration(1000) | StopWhenGradientNormLess(1e-12) | StopWhenStepsizeLess(1e-12)

@time IIF.solveGraphParametric!(fg; stopping_criterion, debug, is_sparse=false, damping_term_min=1e-12, expect_zero_residual=false);

vnd = getVariableSolverData(fg, :b, :parametric)
@test isapprox(vnd.val[1], [0.01, -0.02, 0.03, -0.01, 0.02, -0.03], atol=1e-3)


## Test Case 3

df = CSV.read(@__DIR__() * "/testdata/IMUCase3.csv", DataFrame)

σ_a = 0.037/60*sqrt(100) # noise density m/s²/√Hz - 0.037 m/s/√hr

σ_ω = 0.15*pi*60/180/3600*sqrt(100)  # noise density rad/√Hz - 0.15 deg/√hr
Σy  = diagm([ones(3)*σ_a^2; ones(3)*σ_ω^2])

a_b = SA[0.,0,0]
ω_b = SA[0.,0,0]

accels = map(SVector, zip(df.nax, df.nay, df.naz))
gyros = map(SVector, zip(df.nwx, df.nwy, df.nwz))

# accels = map(SVector, zip(df.ax, df.ay, df.az))
# gyros = map(SVector, zip(df.wx, df.wy, df.wz))

timestamps = df.time


pidx = range(1; step=100, stop=length(timestamps))
factors = map(zip(pidx[1:end-1], pidx[2:end])) do (fr,to)
    r = range(fr,to-1)
    RoME.IMUDeltaFactor(
        SVector{3,Float64}.(accels[r]),
        SVector{3,Float64}.(gyros[r]),
        ones(length(r))*Ts,
        Σy,
        a_b,
        ω_b
    )
end

times = timestamps[pidx]

fg = initfg()
fg.solverParams.graphinit = false

foreach(enumerate(Nanosecond.(times * 10^9))) do (i,nanosecondtime)
    addVariable!(fg, Symbol("x",i-1), Pose3; nanosecondtime)
    addVariable!(fg, Symbol("v",i-1), Position3; nanosecondtime)
end

addFactor!(fg, [:x0], PriorPose3(MvNormal([1.0, 1, 1, 0, 0, 0], diagm(ones(6)*1e-3^2))))
addFactor!(fg, [sortDFG(ls(fg, r"^x"))[end]], PriorPose3(MvNormal([1.0, 1, 1.0, 0, 0, 0], diagm(ones(6)*1e-3^2))))

for l in ls(fg, r"v\d")[1:1]
    addFactor!(fg, [l], PriorPoint3(MvNormal([0.0, 0, 0], diagm(ones(3)*1e-3))))
end

# Variable for bias 
addVariable!(fg, :b, Position{6})

for (i,fac) in enumerate(factors)
    frx = Symbol("x",i-1)
    frv = Symbol("v",i-1)
    tox = Symbol("x",i)
    tov = Symbol("v",i)
    if exists(fg, :b)
        addFactor!(fg, [frx,frv, tox, tov, :b], fac)
    else
        addFactor!(fg, [frx,frv, tox, tov], fac)
    end
end


debug = [:Iteration, :Change,  " | ", :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :damping_term, "\n", :Stop]
debug = [:Iteration, :Cost, " | ", :Stepsize," | ", :GradientNorm," | ", :last_step_successful, "\n", :Stop]
stopping_criterion=StopAfterIteration(1000) | StopWhenGradientNormLess(1e-12) | StopWhenStepsizeLess(1e-12)

@time IIF.solveGraphParametric!(fg; stopping_criterion, debug, is_sparse=false, damping_term_min=1e-12, expect_zero_residual=false);

vnd = getVariableSolverData(fg, :b, :parametric)

@test isapprox(vnd.val[1], [0.01, -0.02, 0.03, -0.01, 0.02, -0.03], atol=1e-3)
