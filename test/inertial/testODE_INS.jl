# Basic test of DERelative


using Test
using DifferentialEquations
using IncrementalInference
using Dates
using Statistics
using TensorCast

## plotting functions

# using Plots
# using Cairo, RoMEPlotting
# Gadfly.set_default_plot_size(25cm,20cm)

##

# @testset "First order DERelative" begin

##

# a user specified ODE in standard form
# inplace `xdot = f(x, u, t)`
# if linear, `xdot = F*x(t) + G*u(t)`
function insKinematic!(dstate, state, u, t)
  # x is a VelPose3 point (assumed ArrayPartition)
  # u is IMU input (assumed [rate; accel])
  Mt = TranslationGroup(3)
  Mr = SpecialOrthogonal(3)
  # convention
  # b is real time body
  # b1 is one discete timestep in the past, equivalent to `r_Sk = r_Sk1 + r_dSk := r_S{k-1} + r_dSk`


  # Using robotics frame fwd-std-dwn <==> North-East-Down
  # ODE cross check taken from Farrell 2008, section 11.2.1, p.388
  # NOTE, Farrell 2008 has gravity logic flipped in some of his equations.
  # WE TAKE gravity as up is positive (think a pendulum hanging back during acceleration)
  # expected gravity in FSD frame (akin to NED).  This is a model of gravity we expect to measure.
  i_G = [0; 0; -9.81] 

  # attitude computer
  w_R_b = state.x[2] # Rotation element
  i_R_b = w_R_b
  # assume body-frame := imu-frame
  b_Ωbi = hat(Mr, Identity(Mr), u[].gyro(t)) # so(3): skew symmetric Lie algebra element
  # assume perfect measurement, i.e. `i` here means measured against native inertial (no coriolis, transport rate, error corrections)
  i_Ṙ_b = i_R_b * b_Ωbi
  # assume world-frame := inertial-frame
  w_Ṙ_b = i_Ṙ_b
  # tie back to the ODE solver

  dstate.x[2] .= w_Ṙ_b
  # Note closed form post integration result (remember exp is solution to first order diff equation)
  # w_R_b = exp(Mr, w_R_b1, b_Ωbi)
  
  # measured inertial acceleration
  b_Abi = u[].accel(t) # already a tangent vector
  # inertial (i.e. world) velocity-dot (accel) by compensating (i.e. removing) expected gravity measurement
  i_V̇ = i_R_b * b_Abi - i_G
  # assume world is inertial frame
  w_V̇ = i_V̇
  dstate.x[3] .= w_V̇ # velocity state
  
  # position-dot (velocity)
  w_V = state.x[3]
  i_V = w_V
  i_Ṗ = i_V
  w_Ṗ = i_Ṗ
  dstate.x[1] .= w_Ṗ
  
  # TODO add biases, see RoME.InertialPose
  # state[4] := gyro bias
  # state[5] := acce bias

  nothing
end

dt = 0.01
N = 101
gyros = [[0.01, 0.0, 0.0] for _ = 1:N]
a0 = [0, 0, -9.81]
accels = [a0]
w_R_b = [1. 0 0; 0 1 0; 0 0 1]
M = SpecialOrthogonal(3)
for g in gyros[1:end-1]
  X = hat(M, Identity(M), g)
  exp!(M, w_R_b, w_R_b, X*dt)
  push!(accels, w_R_b' * a0)
end

gyros_t = linear_interpolation(range(0; step=dt, length=N), gyros)
accels_t = linear_interpolation(range(0; step=dt, length=N), accels)

p = (gyro=gyros_t, accel=accels_t)

p.accel(0.9)

u0 = ArrayPartition([0.0,0,0], Matrix(getPointIdentity(SpecialOrthogonal(3))), [0.,0,0])
tspan = (0.0, 1.0)

prob = ODEProblem(insKinematic!, u0, tspan, Ref(p))

sol = solve(prob)
last(sol)


@test isapprox(last(sol).x[1], [0,0,0]; atol=0.001)
@test isapprox(M, last(sol).x[2], w_R_b; atol=0.001)
@test isapprox(last(sol).x[3], [0,0,0]; atol=0.001)


##


##

# function insTangentFrame!(dstate, state, u, t)
#   Mt = TranslationGroup(3)
#   Mr = SpecialOrthogonal(3)
  
#   #FIXME
#   t_v = [0.,0,0]
#   b_Ωie =  hat(Mr, Identity(Mr), [0.,0,0])

#   t_g = [0; 0; -9.81] 

#   t_R_b = state.x[2] # Rotation element
#   b_Ωib = hat(Mr, Identity(Mr), u[].gyro(t)) # so(3): skew symmetric Lie algebra element
#   t_Ṙ_b = t_R_b * (b_Ωib - b_Ωie)
#   dstate.x[2] .= t_Ṙ_b
#   b_Aib = u[].accel(t) # already a tangent vector
#   t_V̇e = t_R_b * b_Aib - t_g - 2*b_Ωie*t_v
#   dstate.x[3] .= t_V̇e # accel state
#   t_Ve = state.x[3]
#   t_Ṗ = t_Ve
#   dstate.x[1] .= t_Ṗ
#   nothing
# end


@error("WIP BELOW")

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()
# the starting points and "0 seconds"
# `accurate_time = trunc(getDatetime(var), Second) + (1e-9*getNstime(var) % 1)`
addVariable!(fg, :x0, VelPose3, timestamp=DateTime(2000,1,1,0,0,0)) 
# pin with a simple prior
addFactor!(fg, [:x0], Prior(Normal(1,0.01)))

doautoinit!(fg, :x0)

prev = :x0

for i in 1:3

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, Position{1}, timestamp=DateTime(2000,1,1,0,0,5*i))
  # build factor against manifold Manifolds.TranslationGroup(1)
  ode_fac = IIF.DERelative(fg, [prev; nextSym], 
                        Position{1}, 
                        firstOrder!,
                        tstForce,
                        dt=0.05, 
                        problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], ode_fac, graphinit=false )
  initVariable!(fg, nextSym, [zeros(1) for _ in 1:100])

  prev = nextSym
end


## basic sample test

meas = sampleFactor(fg, :x0x1f1, 10)
@test size(meas[1][1],1) == 1
@test size(meas,1) == 10


## do all forward solutions

pts = sampleFactor(fg, :x0f1, 100)

initVariable!(fg, :x0, pts)
pts_ = approxConv(fg, :x0x1f1, :x1)
@cast pts[i,j] := pts_[j][i]
@test 0.3 < Statistics.mean(pts) < 0.4


## check that the reverse solve also works

initVariable!(fg, :x1, pts_)
pts_ = approxConv(fg, :x0x1f1, :x0)
@cast pts[i,j] := pts_[j][i]

# check the reverse solve to be relatively accurate
ref_ = (getBelief(fg, :x0) |> getPoints)
@cast ref[i,j] := ref_[j][i]
@test (pts - ref) |> norm < 1e-4


##

oder_ = DERelative( fg, [:x0; :x3], 
                    Position{1}, 
                    firstOrder!,
                    tstForce, 
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0]
sl = DifferentialEquations.solve(oder_.forwardProblem)

##


# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",layout=(1,1))

# for lb in [:x0; :x1;:x2;:x3]
#   x = getTimestamp(getVariable(fg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [0;1]
#   Plots.plot!(xx, yy, show=true)
# end


##


tfg = initfg()
pts_ = approxConv(fg, :x0f1, :x3, setPPE=true, tfg=tfg)
# initVariable!(tfg, :x3, pts)


##

@cast pts[i,j] := pts_[j][i]

@test getPPE(tfg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(tfg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(tfg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test       Statistics.mean(pts) - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix)[1] < 1.0


##

# plotKDE(tfg, [:x0;:x1;:x2;:x3])


## Now test a full solve

solveTree!(fg);


##


@test getPPE(fg, :x0).suggested - sl(getVariable(fg, :x0) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(fg, :x1).suggested - sl(getVariable(fg, :x1) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(fg, :x2).suggested - sl(getVariable(fg, :x2) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1
@test getPPE(fg, :x3).suggested - sl(getVariable(fg, :x3) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.1


##

end


##

@testset "Damped Oscillator DERelative" begin

## setup some example dynamics

# Lets build an damped oscillator to demonstrate the process in state space
# https://en.wikipedia.org/wiki/Harmonic_oscillator
# ddx/ddt = β dx/dt  -  ω x  +  force[t]
# dx/dt   = dx/dt
function dampedOscillator!(dstate, state, force, t)
  ω = 0.7
  β = -0.3
  dstate[2] = β*state[2] - ω*state[1] + force(t)
  dstate[1] = state[2]
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()

# the starting points and "0 seconds"
addVariable!(fg, :x0, Position{2}, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(MvNormal([1;0],0.01*diagm(ones(2)))))



##

prev = :x0
DT = 2

for i in 1:7

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, Position{2}, timestamp=DateTime(2000,1,1,0,0,DT*i))
  oder = DERelative( fg, [prev; nextSym], 
                      Position{2}, 
                      dampedOscillator!,
                      tstForce, 
                      # (state, var)->(state[1] = var[1]),
                      # (var, state)->(var[1] = state[1]),
                      dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev;nextSym], oder )

  prev = nextSym
end


## check forward and backward solving

pts_ = approxConv(fg, :x0f1, :x0)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [1;0]) < 0.3

initVariable!(fg, :x0, pts_)
X0_ = deepcopy(pts)

pts_ = approxConv(fg, :x0x1f1, :x1)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [0;-0.6]) < 0.4

# now check the reverse direction solving
initVariable!(fg, :x1, pts_)
pts_ = approxConv(fg, :x0x1f1, :x0)
@cast pts[i,j] := pts_[j][i]

@test (X0_ - pts) |> norm < 1e-4


##

tfg = initfg()
for s in ls(fg)
  initVariable!(fg, s, [zeros(2) for _ in 1:100])
end

pts = approxConv(fg, :x0f1, :x7, setPPE=true, tfg=tfg)
# initVariable!(tfg, :x7, pts)



##

# plotKDE(tfg, ls(fg) |> sortDFG, dims=[1] )

##


oder_ = DERelative( fg, [:x0; :x7], 
                    Position{2}, 
                    dampedOscillator!,
                    tstForce, 
                    # (state, var)->(state[1] = var[1]),
                    # (var, state)->(var[1] = state[1]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0;0.0]
sl = DifferentialEquations.solve(oder_.forwardProblem)


## check the solve values are correct


for sym = ls(tfg)
  @test getPPE(tfg, sym).suggested - sl(getVariable(fg, sym) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.2
end


##



# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

# for lb in sortDFG(ls(fg))
#   x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [-1;1]
#   Plots.plot!(xx, yy, show=true)
# end


##

@error "Disabling useMsgLikelihood for DERelative test, follow fix on #1010 as rough guide"
getSolverParams(fg).useMsgLikelihoods = false

solveTree!(fg);


## 


for sym = ls(fg)
  @test getPPE(fg, sym).suggested - sl(getVariable(fg, sym) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.2
end


##

end





##

@testset "Parameterized Damped Oscillator DERelative" begin

## setup some example dynamics

# Lets build an damped oscillator to demonstrate the process in state space
# https://en.wikipedia.org/wiki/Harmonic_oscillator
# ddx/ddt = β dx/dt  -  ω x  +  force[t]
# dx/dt   = dx/dt
# force_ωβ = (data, ωβ)
function dampedOscillatorParametrized!(dstate, state, force_ωβ, t)
  # 3rd variable in this factor graph test example
  force = force_ωβ[1]
  ω     = force_ωβ[2][1]
  β     = force_ωβ[2][2]
  # classic ODE between first and second fg variables
  dstate[2] = β*state[2] - ω*state[1] + force(t)
  dstate[1] = state[2]
  nothing
end

# testing function parameter version (could also be array of data)
tstForce(t) = 0


## build a representative factor graph with ODE built inside

fg = initfg()

# the starting points and "0 seconds"
addVariable!(fg, :x0, Position{2}, timestamp=DateTime(2000,1,1,0,0,0))
# pin with a simple prior
addFactor!(fg, [:x0], Prior(MvNormal([1;0],0.01*diagm(ones(2)))))
doautoinit!(fg, :x0)

# and the new parameterized variable
ω = 0.7
β = -0.3

# these are the stochastic parameters
addVariable!(fg, :ωβ, Position{2}) # timestamp should not matter
# pin with a simple prior
addFactor!(fg, [:ωβ], Prior(MvNormal([ω;β],0.0001*diagm(ones(2)))))
doautoinit!(fg, :ωβ)


##

prev = :x0
DT = 2

for i in 1:7

  nextSym = Symbol("x$i")

  # another point in the trajectory 5 seconds later
  addVariable!(fg, nextSym, Position{2}, timestamp=DateTime(2000,1,1,0,0,DT*i))
  oder = DERelative( fg, [prev; nextSym; :ωβ], 
                      Position{2}, 
                      dampedOscillatorParametrized!,
                      tstForce, # this is passed in as `force_ωβ[1]`
                      # (state, var)->(state[1] = var[1]),
                      # (var, state)->(var[1] = state[1]),
                      # dt=0.05, 
                      problemType=ODEProblem )
  #
  addFactor!( fg, [prev; nextSym; :ωβ], oder, graphinit=false, inflation=0.01 )

  prev = nextSym
end


## check forward and backward solving

pts_ = approxConv(fg, :x0f1, :x0)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [1;0]) < 0.3

initVariable!(fg, :x0, pts_)
X0_ = deepcopy(pts)

pts_ = approxConv(fg, :x0x1ωβf1, :x1)
@cast pts[i,j] := pts_[j][i]
@test norm(Statistics.mean(pts, dims=2) - [0;-0.6]) < 0.4

# now check the reverse direction solving
initVariable!(fg, :x1, pts_)

# failing here
pts_ = approxConv(fg, :x0x1ωβf1, :x0)
@cast pts[i,j] := pts_[j][i]

@test (X0_ - pts) |> norm < 1e-2


##

tfg = initfg()
for s in ls(fg)
  initVariable!(fg, s, [zeros(2) for _ in 1:100])
end

# must initialize the parameters
pts = approxConv(fg, :ωβf1, :ωβ)
initVariable!(fg, :ωβ, pts)

# project forward
forcepath = [:x0f1;]
push!(forcepath, :x0) 
push!(forcepath, :x0x1ωβf1) 
push!(forcepath, :x1)
push!(forcepath, :x1x2ωβf1)
push!(forcepath, :x2)
push!(forcepath, :x2x3ωβf1)
push!(forcepath, :x3)
push!(forcepath, :x3x4ωβf1)
push!(forcepath, :x4)
push!(forcepath, :x4x5ωβf1)
push!(forcepath, :x5)
push!(forcepath, :x5x6ωβf1)
push!(forcepath, :x6)
push!(forcepath, :x6x7ωβf1)
push!(forcepath, :x7)
pts = approxConv(fg, :x0f1, :x7, setPPE=true, tfg=tfg, path=forcepath)


##

# plotKDE(tfg, ls(tfg) |> sortDFG, dims=[1] )


##

# getBelief(fg, :ωβ) |> getPoints

# plotKDE(tfg, :ωβ)

##


oder_ = DERelative( fg, [:x0; :x7; :ωβ], 
                    Position{2}, 
                    dampedOscillatorParametrized!,
                    tstForce,
                    # (state, var)->(state[1] = var[1]),
                    # (var, state)->(var[1] = state[1]),
                    dt=0.05, 
                    problemType=ODEProblem )

oder_.forwardProblem.u0 .= [1.0;0.0]
oder_.data[2] .= [ω;β]
sl = DifferentialEquations.solve(oder_.forwardProblem)



## check the approxConv is working right


for sym in setdiff(ls(tfg), [:ωβ])
  @test getPPE(tfg, sym).suggested - sl(getVariable(fg, sym) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.2
end


## 


# Plots.plot(sl,linewidth=2,xaxis="unixtime [s]",label=["ω [rad/s]" "θ [rad]"],layout=(2,1))

# for lb in sortDFG(ls(fg))
#   x = getTimestamp(getVariable(tfg, lb)) |> DateTime |> datetime2unix
#   xx = [x;x]
#   yy = [-1;1]
#   Plots.plot!(xx, yy, show=true)
# end


## test convolution to the parameter (third) variable

# easy test with good starting points
pts = approxConv(fg, :ωβf1, :ωβ)
initVariable!(fg, :ωβ, pts)

# make sure the other variables are in the right place
pts_ = getBelief(fg, :x0) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [1;0] |> norm < 0.1
pts_ = getBelief(fg, :x1) |> getPoints
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [0;-0.6] |> norm < 0.2


pts_ = approxConv(fg, :x0x1ωβf1, :ωβ)
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [0.7;-0.3] |> norm < 0.1

##

# repeat with more difficult starting point

initVariable!(fg, :ωβ, [zeros(2) for _ in 1:100])

pts_ = approxConv(fg, :x0x1ωβf1, :ωβ)
@cast pts[i,j] := pts_[j][i]
@test Statistics.mean(pts, dims=2) - [0.7;-0.3] |> norm < 0.1


@warn "n-ary DERelative test on :ωβ requires issue #1010 to be resolved first before being reintroduced."
# ## do a complete solve (must first resolve #1010)

# solveTree!(fg);

# ## Solve quality might not yet be good enough for this particular test case

# @test getPPE(fg, :ωβ).suggested - [0.7;-0.3] |> norm < 0.2

# for sym in setdiff(ls(tfg), [:ωβ])
#   @test getPPE(fg, sym).suggested - sl(getVariable(fg, sym) |> getTimestamp |> DateTime |> datetime2unix) |> norm < 0.2
# end


##

end





@error "DERelative not tested for `multihypo=` case yet, see issue #1025"




#