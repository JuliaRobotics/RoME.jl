## Exports


export AbstractSLAM
export ManageSolveSettings, HandlerStateMachine, SolverStateMachine
export manageSolveTree!, stopManageSolveTree!, blockSolvingInProgress
export SLAMCommonHelper, SLAMWrapperLocal
export checkSolveStrideTrigger!, checkSolveStride, triggerSolve!, blockProgress


## types


abstract type AbstractSLAM end

@enum HandlerStateMachine HSMReady HSMHandling HSMOverlapHandling HSMBlocking
@enum SolverStateMachine SSMReady SSMConsumingSolvables SSMSolving


"""
$(TYPEDEF)

DevNotes
- TODO consolidate with SLAMWrapperLocal
"""
mutable struct SLAMWrapper <: AbstractSLAM
  fg::AbstractDFG
  tree::NothingUnion{IncrementalInference.AbstractBayesTree}
  lndmidx::Int
  lastposesym::Symbol
  lastpose::SE3
  lbl2usrid::Dict{Symbol,Int}
  usrid2lbl::Dict{Int, Symbol}
  SLAMWrapper() = new()
  SLAMWrapper(a,b,c) = new(a,b,c,:x1, SE3(0), Dict{Symbol,Int}(), Dict{Int,Symbol}()) # TODO to be deprecated
end

"""
    $TYPEDEF

Container for settings and synchronization tools for use with manageSolveTree.
"""
mutable struct ManageSolveSettings
  solveStride
  loopSolver
  solvables
  solveInProgress
  poseSolveToken
  canTakePoses::Condition
  drtCurrent
  # constructor
  ManageSolveSettings(;solveStride=10,
                       loopSolver=true,
                       solvables=Channel{Vector{Symbol}}(100),
                       solveInProgress=SSMReady,
                       poseSolveToken=Channel{Int}(2), #  ensure only one solve can occur at a time (see blockProgress)
                       canTakePoses::Condition=Condition(),
                       drtCurrent=(:null,:null)) = new(solveStride, loopSolver, solvables, solveInProgress, poseSolveToken, canTakePoses, drtCurrent)
end


mutable struct SLAMCommonHelper
  lastPoseOdomBuffer::SE3 # common helper
  SLAMCommonHelper(lpo::SE3=SE3(0)) = new(lpo)
end


mutable struct SLAMWrapperLocal{G <: AbstractDFG} <: AbstractSLAM
  dfg::G
  poseCount::Int
  frameCounter::Int
  poseStride::Int # pose every frameStride (naive pose trigger)
  helpers::SLAMCommonHelper
  solveSettings::ManageSolveSettings
end


SLAMWrapperLocal(;dfg::G=initfg(),
                  poseCount::Int=0,
                  frameCounter::Int=0,
                  poseStride::Int=10,
                  helpers::SLAMCommonHelper=SLAMCommonHelper(),
                  solveSettings::ManageSolveSettings=ManageSolveSettings() ) where {G <: AbstractDFG}= SLAMWrapperLocal{G}(dfg,
                      poseCount, frameCounter, poseStride, helpers, solveSettings)
#

"""
    $SIGNATURES

Trigger a factor graph `solveTree!(slam.dfg,...)` after clearing the solvable buffer `slam.??` (assuming the `manageSolveTree!` task is already running).

Notes
- Used in combination with `manageSolveTree!`
"""
triggerSolve!(slam::SLAMWrapperLocal) = put!(slam.solveSettings.poseSolveToken, slam.poseCount)

function checkSolveStride(slam::SLAMWrapperLocal)
  slam.poseCount % slam.solveSettings.solveStride == 0
end

"""
    $SIGNATURES

Check and trigger a solve if `slam.poseCount` reached the solve stride `slam.solveSettings.solveStride`.

Notes
- Used in combination with `manageSolveTree!`

Related

triggerSolve!
"""
function checkSolveStrideTrigger!(slam::SLAMWrapperLocal; force::Bool=false)
  if force || checkSolveStride(slam)
    @info "trigger a solve, $(length(slam.solveSettings.poseSolveToken.data)) ====================================="
    triggerSolve!(slam)
    @info "after slam solve trigger"
    # also update DRT containers

    return true
  end
  return false
end


"""
    $SIGNATURES

Block progress on current task until SLAMWrapper is ready to continue taking on more data.

Notes:
- Use to block addition of more data into slam.dfg object, usually used for preventing duplicate solutions being invoked concurrently.
  - Consider dfg segments: 1. already busy solving, 2. recently added; and now attempting to launch solve on 1&2 before solve on 1. finishes.
- Also see dead reckon tether (DRT).
- Assumes unsolved poses are always added and never removed before being solved.

Related

manageSolveTree!, SLAMWrapperLocal
"""
function blockProgress(slam::SLAMWrapperLocal)
  mss = slam.solveSettings
  # determine if solve is in progress
  # mss.solveInProgress

  # isready means solve is in progress
  if 1 < length(mss.poseSolveToken.data) && checkSolveStride(slam)
    @warn "blockProgress called and forced to wait on previous slam solve not having completed."
    wait(mss.canTakePoses)
  end
end

"""
    $SIGNATURES

Stops a `manageSolveTree!` session.  Usually up to the user to do so as a SLAM process comes to completion.

Related

manageSolveTree!
"""
function stopManageSolveTree!(slam::SLAMWrapperLocal)
  slam.solveSettings.loopSolver = false
end

"""
    $SIGNATURES

Block the progression of calling task if `::SLAMWrapperLocal` is being solved by a presumed `manageSolveTree!` task.
"""
function blockSolvingInProgress(slam::SLAMWrapperLocal)
  if slam.solveSettings.solveInProgress != SSMReady
    @info "blockSolvingInProgress on ::SLAMWrapperLocal"
    wait(slam.solveSettings.canTakePoses)
    @info "blockSolvingInProgress notified of completion."
  end
end

"""
    $SIGNATURES

Asynchronous solver manager that can run concurrently while other Tasks are modifying a common distributed factor graph object.

Notes
- When adding Variables and Factors, use `solvable=0` to disable the new fragments until ready for inference.
  - e.g. `addVariable!(fg, :x45, Pose2, solvable=0)`
  - These parts of the factor graph can simply be activated for solving `setSolvable!(fg, :x45, 1)`
"""
function manageSolveTree!(dfg::AbstractDFG,
                          mss::ManageSolveSettings;
                          dbg::Bool=false,
                          timinglog=Base.stdout,
                          limitfixeddown::Bool=true  )
  #
  @info "logpath=$(getLogPath(dfg))"

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
    while mss.loopSolver
      t0 = time_ns()
      solvecycle += 1
      # add any newly solvables (atomic)
      while !isready(mss.solvables) && mss.loopSolver && !isready(mss.poseSolveToken)
        @info "waiting on !isready(mss.solvables)=$(!isready(mss.solvables)) and mss.loopSolver=$(mss.loopSolver)"
        sleep(0.2)
      end
      dt_wait = (time_ns()-t0)/1e9

      #add any new solvables
      while isready(mss.solvables) && mss.loopSolver
        mss.solveInProgress = SSMConsumingSolvables
        @show tosolv = take!(mss.solvables)
        for sy in tosolv
          # setSolvable!(dfg, sy, 1) # see DFG #221
          # TODO temporary workaround
          getfnc = occursin(r"f", string(sy)) ? getFactor : getVariable
          getfnc(dfg, sy).solvable = 1
        end
      end
      dt_solvable = (time_ns()-t0)/1e9

      @info "Ensure all new variables initialized"
      initAll!(dfg)
      dt_init = (time_ns()-t0)/1e9
      dt_disengage = 0.0
      dt_save1 = 0.0
      dt_solve = 0.0

      # mss.solveInProgress = SSMReady

      # solve only every 10th pose
      @show length(mss.poseSolveToken.data)
      if 0 < length(mss.poseSolveToken.data)
      # if 10 <= mss[:poseStride]
        @info "reduce problem size by disengaging older parts of factor graph"
        setSolvableOldPoses!(dfg, youngest=getSolverParams(dfg).qfl+round(Int,mss.solveStride/2), oldest=100, solvable=0)
        dt_disengage = (time_ns()-t0)/1e9

        # set up state machine flags to allow overlapping or block
        mss.solveInProgress = SSMSolving # obsolete
        # mss[:poseStride] = 0

        # do the actual solve (with debug saving)
        lasp = getLastPoses(dfg, filterLabel=r"x\d", number=1)[1]
        !dbg ? nothing : saveDFG(dfg, joinpath(getLogPath(dfg), "fg_before_$(lasp)"))
        dt_save1 = (time_ns()-t0)/1e9
        # constrain solve with the latest pose at the top
        # @show latestPose = intersect(getLastPoses(dfg, filterLabel=r"x\d", number=12), ls(dfg, r"x\d", solvable=1))[end]
        tree = solveTree!(dfg, tree) # , variableConstraints=[latestPose;]
        dt_solve = (time_ns()-t0)/1e9
        !dbg ? nothing : saveDFG(dfg, joinpath(getLogPath(dfg), "fg_after_$(lasp)"))

        # unblock LCMLog reader for next STRIDE segment
        mss.solveInProgress = SSMReady

        # adjust latest RTT after solve, latest solved -- hard coded pose stride 10
        lastList = sortDFG(ls(dfg, r"x\d+9\b|x9\b", solvable=1))
        if 0 < length(lastList)
          lastSolved = lastList[end]
          mss.drtCurrent = (lastSolved, Symbol("drt_"*string(lastSolved)[2:end]))
        end

        # remove a token to allow progress to continue
        gotToken = take!(mss.poseSolveToken)

        # notify any processes that might be waiting on current solve to complete
        # lock(mss.canTakePoses)  # uncomment when upgrading to Threads.Condition
        notify(mss.canTakePoses) # = HSMHandling
        # unlock(mss.canTakePoses)

        "end of solve cycle, token=$gotToken" |> println
      else
        "sleep a solve cycle" |> println
        mss.solveInProgress = SSMReady
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
