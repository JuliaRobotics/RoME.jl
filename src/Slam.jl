
## types

abstract type AbstractSLAM end

@enum HandlerStateMachine HSMReady HSMHandling HSMOverlapHandling HSMBlocking
@enum SolverStateMachine SSMReady SSMConsumingSolvables SSMSolving


"""
$(TYPEDEF)
"""
mutable struct SLAMWrapper <: AbstractSLAM
  fg::IncrementalInference.FactorGraph
  tree::NothingUnion{IncrementalInference.BayesTree}
  lndmidx::Int
  lastposesym::Symbol
  lastpose::SE3
  lbl2usrid::Dict{Symbol,Int}
  usrid2lbl::Dict{Int, Symbol}
  SLAMWrapper() = new()
  SLAMWrapper(a,b,c) = new(a,b,c,:x1, SE3(0), Dict{Symbol,Int}(), Dict{Int,Symbol}()) # TODO to be deprecated
end




function addOdoFG!(slaml::SLAMWrapper, odo::Pose3Pose3;
                  N::Int=100, solvable::Int=1,
                  saveusrid::Int=-1)
  #
  vprev = getVert(slaml.fg, slaml.lastposesym)
  # vprev, X, nextn = getLastPose(fgl)
  npnum = parse(Int,string(slaml.lastposesym)[2:end]) + 1
  nextn = Symbol("x$(npnum)")
  vnext = addVariable!(slaml.fg, nextn, Pose2(labels=["POSE";]), N=N, solvable=solvable)
  # vnext = addVariable!(slaml.fg, nextn, getVal(vprev) ⊕ odo, N=N, solvable=solvable, labels=["POSE"])
  slaml.lastposesym = nextn
  fact = addFactor!(slaml.fg, [vprev;vnext], odo)

  if saveusrid > -1
    slaml.lbl2usrid[nextn] = saveusrid
    slaml.usrid2lbl[saveusrid] = nextn
  end
  return vnext, fact
end


function addposeFG!(slaml::SLAMWrapper,
      constrs::Vector{IncrementalInference.FunctorInferenceType};
      N::Int=100,
      solvable::Int=1,
      saveusrid::Int=-1   )
  #
  vprev = getVert(slaml.fg, slaml.lastposesym)

  npnum = parse(Int,string(slaml.lastposesym)[2:end]) + 1
  nextn = Symbol("x$(npnum)")
  # preinit
  vnext = nothing
  if !haskey(slaml.fg.IDs, nextn)
    vnext = addVariable!(slaml.fg, nextn, Pose2, N=N, solvable=solvable)
    # vnext = addVariable!(slaml.fg, nextn, getVal(vprev), N=N, solvable=solvable, labels=["POSE"])
  else
    vnext = getVert(slaml.fg, nextn) #, api=localapi # as optimization, assuming we already have latest vnest in slaml.fg
  end
  slaml.lastposesym = nextn

  addsubtype(fgl::FactorGraph, vprev, vnext, cc::IncrementalInference.FunctorPairwise) = addFactor!(fgl, [vprev;vnext], cc)
  addsubtype(fgl::FactorGraph, vprev, vnext, cc::IncrementalInference.FunctorSingleton) = addFactor!(fgl, [vnext], cc)

  facts = Graphs.ExVertex[]
  PP = BallTreeDensity[]
  for cns in constrs
    fa = addsubtype(slaml.fg, vprev, vnext, cns)
    push!(facts, fa)
  end

  # set node val from new constraints as init
  val, = predictbelief(slaml.fg, vnext, facts, N=N)
  setVal!(vnext, val)
  IncrementalInference.dlapi.updatevertex!(slaml.fg, vnext)

  if saveusrid > -1
    slaml.lbl2usrid[nextn] = saveusrid
    slaml.usrid2lbl[saveusrid] = nextn
  end
  return vnext, facts
end



# Requires
#  dashboard: SOLVESTRIDE, loopSolver, solvables, loopSolver, solveInProgress, poseSolveToken, canTakePoses, drtCurrent
#
function manageSolveTree!(dfg::AbstractDFG,
                          dashboard::Dict;
                          dbg::Bool=false,
                          timinglog=Base.stdout,
                          limitfixeddown::Bool=true  )
  #
  @info "logpath=$(getLogPath(dfg))"
  getSolverParams(dfg).drawtree = true
  getSolverParams(dfg).qfl = dashboard[:SOLVESTRIDE]
  getSolverParams(dfg).isfixedlag = true
  getSolverParams(dfg).limitfixeddown = limitfixeddown

  # allow async process
  # getSolverParams(dfg).async = true

  # prep with empty tree
  tree = emptyBayesTree()

  # needs to run asynchronously
  ST = @async begin
    while @show length(ls(dfg, :x0, solvable=1)) == 0
      "waiting for prior on x0" |> println
      sleep(1)
    end
    solvecycle = 0
    # keep solving
    while dashboard[:loopSolver]
      t0 = time_ns()
      solvecycle += 1
      # add any newly solvables (atomic)
      while !isready(dashboard[:solvables]) && dashboard[:loopSolver]
        sleep(0.2)
      end
      dt_wait = (time_ns()-t0)/1e9

      #add any new solvables
      while isready(dashboard[:solvables]) && dashboard[:loopSolver]
        dashboard[:solveInProgress] = SSMConsumingSolvables
        @show tosolv = take!(dashboard[:solvables])
        for sy in tosolv
          # setSolvable!(dfg, sy, 1) # see DFG #221
          # TODO temporary workaround
          getfnc = occursin(r"f", string(sy)) ? getFactor : getVariable
          getfnc(dfg, sy).solvable = 1
        end
      end
      dt_solvable = (time_ns()-t0)/1e9

      @info "Ensure all new variables initialized"
      ensureAllInitialized!(dfg)
      dt_init = (time_ns()-t0)/1e9
      dt_disengage = 0.0
      dt_save1 = 0.0
      dt_solve = 0.0

      dashboard[:solveInProgress] = SSMReady

      # solve only every 10th pose
      if 0 < length(dashboard[:poseSolveToken].data)
      # if 10 <= dashboard[:poseStride]
        @info "reduce problem size by disengaging older parts of factor graph"
        setSolvableOldPoses!(dfg, youngest=getSolverParams(dfg).qfl+round(Int,dashboard[:SOLVESTRIDE]/2), oldest=100, solvable=0)
        dt_disengage = (time_ns()-t0)/1e9

        # set up state machine flags to allow overlapping or block
        dashboard[:solveInProgress] = SSMSolving
        # dashboard[:poseStride] = 0

        # do the actual solve (with debug saving)
        lasp = getLastPoses(dfg, filterLabel=r"x\d", number=1)[1]
        !dbg ? nothing : saveDFG(dfg, joinpath(getLogPath(dfg), "fg_before_$(lasp)"))
        dt_save1 = (time_ns()-t0)/1e9
        # constrain solve with the latest pose at the top
        # @show latestPose = intersect(getLastPoses(dfg, filterLabel=r"x\d", number=12), ls(dfg, r"x\d", solvable=1))[end]
        tree, smt, hist = solveTree!(dfg, tree, maxparallel=1000) # , variableConstraints=[latestPose;]
        dt_solve = (time_ns()-t0)/1e9
        !dbg ? nothing : saveDFG(dfg, joinpath(getLogPath(dfg), "fg_after_$(lasp)"))

        # unblock LCMLog reader for next STRIDE segment
        dashboard[:solveInProgress] = SSMReady
        # de-escalate handler state machine
        dashboard[:canTakePoses] = HSMHandling

        # adjust latest RTT after solve, latest solved -- hard coded pose stride 10
        lastList = sortDFG(ls(dfg, r"x\d+9\b|x9\b", solvable=1))
        if 0 < length(lastList)
          lastSolved = lastList[end]
          dashboard[:drtCurrent] = (lastSolved, Symbol("drt_"*string(lastSolved)[2:end]))
        end

        # remove a token to allow progress to continue
        gotToken = take!(dashboard[:poseSolveToken])
        "end of solve cycle, token=$gotToken" |> println
      else
        "sleep a solve cycle" |> println
        sleep(0.2)
      end
      dt_finish = (time_ns() - t0)/1e9

      # store timing results
      println(timinglog, "$solvecycle, $t0, $dt_wait, $dt_solvable, $dt_init, $dt_disengage, $dt_save1, $dt_solve, $(tree.buildTime), $dt_finish")
    end
  end
  return ST
end





#
