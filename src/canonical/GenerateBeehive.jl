# generate a randomized Beehive



"""
    $SIGNATURES

Pretend a bee is walking in a hive where each step (pose) follows one edge of an imaginary honeycomb lattice, 
and at after each step a new direction left or right is stochastically chosen and the process repeats.

Notes
- The keyword `locality=1` is a positive `::Real` ∈ [0,∞) value, where higher numbers imply direction decisions are more sticky for multiple steps.
- Use keyword callback function `postpose_cb = (fg, lastpose) -> ...` to hook in your own features right after each new pose step.

DevNotes
- TODO rewrite as a recursive generator function instead.

See also: [`generateGraph_Honeycomb!`](@ref), [`generateGraph_Hexagonal`](@ref), [`generateGraph_ZeroPose`](@ref)
"""
function generateGraph_Beehive!(poseCountTarget::Int=10;
                                graphinit::Bool = true,
                                dfg::AbstractDFG = LocalDFG{SolverParams}(solverParams=SolverParams(graphinit=graphinit)),
                                useMsgLikelihoods::Bool=getSolverParams(dfg).useMsgLikelihoods,
                                solvable::Int = 1,
                                refKey::Symbol = :simulated,
                                addLandmarks::Bool = true,
                                landmarkSolvable::Int=0,
                                poseRegex::Regex=r"x\d+",
                                pose0::Symbol=Symbol(match(r"[a-zA-Z_]+", poseRegex.pattern).match, 0),
                                yaw0::Real = ([0.0;-2pi/3;2pi/3])[rand(1:3)],
                                μ0::AbstractVector{<:Real} = [0;0;yaw0],                                  
                                postpose_cb::Function=(fg_,latestpose)->(),
                                locality::Real=1,
                                atol::Real=1     )
  #

  # does anything exist in the graph yet
  posecount = if pose0 in ls(dfg)
    # what is the last pose
    lastPose = (ls(dfg, poseRegex) |> sortDFG)[end]
    # get latest posecount number
    match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  else
    # pick a random staring direction
    # initial zero pose
    generateGraph_ZeroPose(dfg=dfg, varType=Pose2, μ0=μ0, graphinit=graphinit, postpose_cb=postpose_cb)
    getSolverParams(dfg).useMsgLikelihoods = useMsgLikelihoods
    
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(dfg, pose0, refKey=refKey, solvable=landmarkSolvable, atol=atol, graphinit=false, _doHoneycomb=false)

    # staring posecount (i.e. :x0)
    0
  end

  # keep adding poses until the target number is reached
  direction = rand(1:2) === 1 ? :left : :right
  _locality = 1/(1.0+locality)
  while posecount < poseCountTarget
    #change or keep direction according to locality keyword
    direction = rand() < _locality ? (direction == :left ? :right : :left) : direction
    posecount = _offsetHexLeg(dfg, posecount, direction=direction, graphinit=graphinit, 
                              landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget, 
                              atol=atol, postpose_cb=postpose_cb, _doHoneycomb=false)
  end

  # NOTE solvable forced for everything at this time
  setSolvable!.(dfg, ls(dfg),  solvable)
  setSolvable!.(dfg, lsf(dfg), solvable)

  return dfg  
end


#
