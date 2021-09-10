# simulate circle multimodal vs outlier example

# using Revise

using ArgParse



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--altmode"
            help = "fraction for alternative mode"
            arg_type = Float64
            default = 0.7
        "--resolve"
            help = "should a repeat solve be performed in the loop"
            action = :store_true
        "--CYCLES"
            help = "How many cycles to drive"
            arg_type = Int
            default = 10
        "--spreadNH"
            help = "Scale adjustment for spreading nullhypo (entropy)"
            arg_type = Float64
            default = 3.0
    end

    return parse_args(s)
end

pargs = parse_commandline()

using Distributed
# addprocs(8) # use -p8 instead

using Random

using RoME, DistributedFactorGraphs
using RoMEPlotting
using KernelDensityEstimatePlotting
@everywhere using RoME

using FileIO
using JLD2
using JSON2
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
CYCLES = pargs["CYCLES"]
LANDMARKS = 10
altFrac = pargs["altmode"]
mainFrac = 1-pargs["altmode"]

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
if pargs["CYCLES"] == 10
  @assert sum(union(sequence...)) == 55
end

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

getSolverParams(fg).spreadNH = pargs["spreadNH"]

ensureAllInitialized!(fg)

## store arguments in results log

argstr = JSON2.write(pargs)

Base.mkpath(getLogPath(fg))
fid = open(joinLogPath(fg, "args.json"), "w")
JSON2.write(fid, pargs)
close(fid)
fid = open(joinLogPath(fg, "..", "results.log"), "a")
Base.write(fid, "$(getLogPath(fg)), RoME/examples/OutlierVsMultimodal.jl, ")
JSON2.write(fid, pargs)
Base.write(fid, "\n")
close(fid)


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

plfl1 = drawPosesLandms(fg, spscale=1.0,title=getLogPath(fg)*"\n$argstr")

## prepare factors to use

BRdistr = Dict{Symbol,Dict{Symbol,Tuple}}()
for i in 1:2*SIZE
  lmid = Symbol(i <= SIZE ? "l$i" : "l$(i-SIZE)_0")
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
tree = solveTree!(fg)


plfl1 = drawPosesLandms(fg, spscale=1.0,title=getLogPath(fg)*"\n$argstr")


# lookup
##  change to alternate location, one at a time


# activate solvable on alternative landmark locations for remainer of test
for l in 1:10
  lmid_0 = Symbol("l$(l)_0")
  setSolvable!(fg, lmid_0, 1)
end

POSEOFFSET=0

## store

@save joinLogPath(fg, "lookup.jld2") lookup

plfl1 = drawPosesLandms(fg, spscale=1.0,title=getLogPath(fg)*"\n$argstr", landmsPPE=:max, contour=true) #, posesPPE=:max)
plfl1 |> PDF(joinLogPath(fg, "plot_x20_before.pdf"), 20cm, 17cm)

##

# Drive CYCLES-1 more loops
i = 2
for i in 2:CYCLES-1

global fg
global lookup
global SIZE
global BRdistr
global pargs
global POSEOFFSET += 10
# drive the new loop without landmark detections (dont solve yet)
fg = generateCanonicalFG_Circle(i*SIZE, fg=fg, kappaOdo=0.1, loopClosure=false, landmark=false, cyclePoses=10)

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
    addFactor!(fg, [opsid;lmid;lmid_0], Pose2Point2BearingRange(ppbr...), multihypo=[1;mainFrac;altFrac])
  end

  # reinit all lm_0 with the latest factor
  if 0 < length(ls(fg, lmid_0))
    pts = approxConv(fg, ls(fg, lmid_0)[1], lmid_0  )
    initManual!(fg,lmid_0,pts)
  end
end

getSolverParams(fg).dbg = true
saveDFG(fg, joinLogPath(fg, "fg_x$(POSEOFFSET+10)before"))
getSolverParams(fg).maxincidence = 1000
tree = solveTree!(fg, recordcliqs=ls(fg))
saveDFG(fg, joinLogPath(fg, "fg_x$(POSEOFFSET+10)_solve"))

plfl1 = drawPosesLandms(fg, spscale=1.0,title=getLogPath(fg)*"\n$argstr", landmsPPE=:max, contour=true)
plfl1 |> PDF(joinLogPath(fg, "plot_x$(POSEOFFSET+10)_solve.pdf"), 20cm, 17cm)

if pargs["resolve"]
  tree = solveTree!(fg, recordcliqs=ls(fg))
  saveDFG(fg, joinLogPath(fg, "fg_x$(POSEOFFSET+10)_resolve"))
  plfl1 = drawPosesLandms(fg, spscale=1.0,title=getLogPath(fg)*"\n$argstr", landmsPPE=:max, contour=true)
  plfl1 |> PDF(joinLogPath(fg, "plot_x$(POSEOFFSET+10)_resolve.pdf"), 20cm, 17cm)
end

end



## test solve

getSolverParams(fg).dbg = true
getSolverParams(fg).maxincidence = 1000
tree = solveTree!(fg, recordcliqs=ls(fg)) #[:l9_0;:x14])
saveDFG(fg, joinLogPath(fg, "fg_x$(POSEOFFSET+10)_final"))


##

plfl1 = drawPosesLandms(fg, spscale=1.0,title=getLogPath(fg)*"\n$argstr", landmsPPE=:max)
plfl1 |> PDF(joinLogPath(fg, "plot_x$(POSEOFFSET+10)_final.pdf"), 20cm, 17cm)


#
