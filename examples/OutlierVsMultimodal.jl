# simulate circle multimodal vs outlier example

# using Revise

using Distributed
using Random

using RoME, DistributedFactorGraphs
using RoMEPlotting
using KernelDensityEstimatePlotting
@everywhere using RoME

using FileIO
using JLD2
using Gadfly
Gadfly.set_default_plot_size(35cm,25cm)

include(joinpath(@__DIR__, "FixedLagCirclePlotting.jl"))

##  Landmark ground truth

# assume all 10 landmarks changing position from rows 1:10 to 11:20, randomly after each circle
L = [0.0 20.0;      # L1
     30.0 30.0;     # L2
     -10.0 20.0;    # L3
     20.0 30.0;     # L4
     -20.0 30.0;    # L5
     -20.0 10.0;    # L6
     20.0 -0.0;     # L7
     10.0 10.0;     # L8
     0.0 40.0;      # L9
     30.0 10.0;     # L10
     10.0 -10.0;    # L1_0
     10.0 40.0;     # L2_0
     -10.0 -10.0;   # L3_0
     -10.0 30.0;    # L4_0
     30.0 -0.0;     # L5_0
     20.0 20.0;     # L6_0
     0.0 10.0;      # L7_0
     30.0 20.0;     # L8_0
     20.0 -10.0;    # L9_0
     -10.0 10.0]    # L10_0

# working from previous save stored right here in the example folder
# @load string(joinpath(@__DIR__, "LassoBearingRange.jld2")) BR

SIZE = 10
CYCLES = 10
LANDMARKS = 10


## divise landmark rotation schedule


movePerCycle = rand(Categorical([0.6;0.4]),CYCLES)
idx = sum(cumsum(movePerCycle) .< CYCLES)
movePerCycle = movePerCycle[1:idx]
push!(movePerCycle, CYCLES-cumsum(movePerCycle)[end])
order = shuffle(1:CYCLES)
sequence = Vector{Vector{Int}}()
total = 1
for n in 1:length(movePerCycle)
  global total
  push!(sequence, order[total:total+movePerCycle[n]-1] )
  total += movePerCycle[n]
end
@assert sum(union(sequence...)) == 55

sequence

## convert sequence to landmark lookups

lookup = Vector{Vector{Symbol}}()
push!(lookup, [Symbol("l$i") for i in 1:10])

for i in 1:9
  lu = deepcopy(lookup[i])
  # before each loop, change location of some landmarks
  seq = i <= length(sequence) ? sequence[i] : Int[]
  for s in seq
    lu[s] = Symbol("l$(s)_0")
  end

  # push new vector into lookup
  push!(lookup, lu)
end

# this is the sequence of landmark transitions
lookup


# BRdistr[:l1]
# BRdistr[:l1_0]


### SETUP first circle----------------------------------------------

# drive first circle
fg = generateCanonicalFG_Circle(SIZE, kappaOdo=0.1, loopClosure=false, landmark=false, cyclePoses=10)


ensureAllInitialized!(fg)

## add initial landmarks

BR = Dict{Symbol,Dict{Symbol, Tuple}}()
# add the landmarks
for i in 1:2*CYCLES
  lmid = Symbol(i <= 10 ? "l$i" : "l$(i-10)_0")
  # lmid = Symbol("l$i")
  addVariable!(fg, lmid, RoME.Point2)
  pts = rand(MvNormal(L[i,:],diagm([0.01;0.01].^2)),100)
  initManual!(fg, lmid, manikde!(pts, RoME.Point2))
  setVariablePosteriorEstimates!(fg,lmid)

  BR[lmid] = Dict{Symbol,Tuple}()
  varNear, varDist = findVariablesNear(fg, L[i,:], r"x\d", number=10)
  for vsym in varNear[varDist .< 25]
    # calculate bearing and range factor
    bear, rang = predictBodyBR(fg, vsym, lmid)
    if abs(bear) < pi/3 && vsym != :x10
      BR[lmid][vsym] = (bear, rang)
    end
  end
end
for i in CYCLES+1:2*CYCLES
  lmid = Symbol(i <= 10 ? "l$i" : "l$(i-10)_0")
  # vsym = Symbol("l$i")
  setSolvable!(fg, lmid, 0)
  # ( x->deleteFactor!(fg, x) ).( ls(fg, vsym) )  # no factors added yet
  # deleteVariable!(fg, lmid)
end
# @save joinpath(@__DIR__, "LassoBearingRange.jld2") BR


## see what is going on

plfl1 = drawPosesLandms(fg, spscale=1.0)

## prepare factors to use

BRdistr = Dict{Symbol,Dict{Symbol,Tuple}}()
for i in 1:2*CYCLES
  lmid = Symbol(i <= CYCLES ? "l$i" : "l$(i-CYCLES)_0")
  BRdistr[lmid] = Dict{Symbol,Tuple}()
  for (vsym, br) in BR[lmid]
    @show bear, rang = br
    BRdistr[lmid][vsym] = (Normal(bear,0.05),Normal(rang,0.3))
  end
end

BRdistr


## add factors to the landmarks on first cycle

@assert ls(fg, :l1) |> length == 0

OFFSET=0
for id in 1:CYCLES
  lmid = Symbol("l$(id)")
  lmid_0 = Symbol("l$(id)_0")
  for (psid, ppbr) in BRdistr[lmid]
    # global OFFSET
    adjpsid = psid # add offset
    addFactor!(fg, [adjpsid;lmid], Pose2Point2BearingRange(ppbr...))
  end
end


##


getSolverParams(fg).drawtree = true
tree, smt, hist = solveTree!(fg)


plfl1 = drawPosesLandms(fg, spscale=1.0)



##  change to alternate location, one at a time


# activate solvable on alternative landmark locations for remainer of test
for l in 1:10
  lmid_0 = Symbol("l$(l)_0")
  setSolvable!(fg, lmid_0, 1)
end

# Drive CYCLES-1 more loops
POSEOFFSET=0
i = 2
# for i in 2:2 #CYCLES-1

global POSEOFFSET += 10
# drive the new loop without landmark detections (dont solve yet)
fg = generateCanonicalFG_Circle(2*SIZE, fg=fg, kappaOdo=0.1, loopClosure=false, landmark=false, cyclePoses=10)

# add modified landmark sighting measurements
for l in 1:LANDMARKS
  lmidBr = lookup[i][l]
  lmid = Symbol("l$(l)")
  lmid_0 = Symbol("l$(l)_0")
  for (psid, ppbr) in BRdistr[lmidBr]
    # apply pose offset
    psnum = parse(Int, string(psid)[2:end]) + POSEOFFSET
    opsid = Symbol("x$psnum")
    # add the few BR factors from this lmid accordingly
    addFactor!(fg, [opsid;lmid;lmid_0], Pose2Point2BearingRange(ppbr...), multihypo=[1;0.3;0.7])
  end
end



# end


# drawGraph(fg)

## dev testing

fg = generateCanonicalFG_Circle(2*SIZE, fg=fg, kappaOdo=0.1, loopClosure=false, landmark=false, cyclePoses=10)


ensureAllInitialized!(fg)
plfl1 = drawPosesLandms(fg, spscale=1.0)

getSolverParams(fg).dbg = true
tree, smt, hist = solveTree!(fg)


BRdistr[:l7]
BRdistr[:l7_0]


plotLocalProduct(fg, :l7_0, levels=8)

plotLocalProduct(fg, :l7, levels=8)



keep_l7_0 = getKDE(fg, :l7_0) |> deepcopy
initManual!(fg, :l7_0, deepcopy(getKDE(fg, :l7)))


ls(fg, :l7_0)

getSolverData(getFactor(fg, :x16l7l7_0f1))


getFactorType(fg, :x0l7f1).bearing
getFactorType(fg, :x10l7l7_0f1)


ls(fg, :x0)
ls(fg, :l7_0)

pts = approxConv(fg, :x16l7l7_0f1, :l7_0)
pts = approxConv(fg, :x16l7l7_0f1, :l7)
pts = approxConv(fg, :x0l7f1, :l7)


plotKDE(manikde!(pts, Point2))


## Something is wrong with :l3

# Juno.@enter approxConv(fg, :x16l7l7_0f1, :l7)

plotLocalProduct(fg, :l3, levels=8)
plotLocalProduct(fg, :l3_0, levels=8)

pts = approxConv(fg, :x16l3l3_0f1, :l3)
pts = approxConv(fg, :x16l3l3_0f1, :l3_0)


drawTree(tree, show=true)
spyCliqMat(tree, :l3)

plotUpMsgsAtCliq(tree, :x5, :l3)
plotUpMsgsAtCliq(tree, :l7, :l3) # problem is in here perhaps


sfg_x7 = buildCliqSubgraph(fg, tree, :x7)


um1 = getUpMsgs(tree,:x2)
um2 = getUpMsgs(tree,:l6)

addMsgFactors!(sfg_x7,um1)
addMsgFactors!(sfg_x7,um2)

drawGraph(sfg_x7)

plotLocalProduct(sfg_x7,:l3)

#
