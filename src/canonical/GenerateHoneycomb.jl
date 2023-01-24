

global _honeycombRecipe = Dict{Symbol,Symbol}(
  :l6   =>:l0,
  :l18  =>:l0,
  :l42  =>:l0,
  :l57  =>:l0,
  :l63  =>:l0,
  :l70  =>:l0,
  :l10  =>:l1,
  :l50  =>:l1,
  :l56  =>:l1,
  :l69  =>:l1,
  :l43  =>:l2,
  :l49  =>:l2,
  :l61  =>:l2,
  :l36  =>:l3,
  :l54  =>:l3,
  :l28  =>:l4,
  :l34  =>:l4,
  :l47  =>:l4,
  :l12  =>:l5,
  :l19  =>:l5,
  :l26  =>:l5,
  :l33  =>:l5,
  :l40  =>:l5,
  :l13  =>:l7,
  :l25  =>:l7,
  :l17  =>:l8,
  :l58  =>:l9,
  :l64  =>:l9,
  :l35  =>:l11,
  :l41  =>:l11,
  :l48  =>:l11,
  :l55  =>:l11,
  :l62  =>:l11,
  :l20  =>:l14,
  :l32  =>:l14,
  :l24  =>:l15,
  :l27  =>:l21,
  :l39  =>:l21,
  :l31  =>:l22,
  :l38  =>:l29,
  :l46  =>:l37,
  :l53  =>:l44,
  :l60  =>:l51,

  # pose offset legs
  :x41  =>:left,
  :x63  =>:left,
  :x78  =>:left,
)

# piecewise posecount construction
# 7,14,21,28,35,43,50



function _addLandmarkBeehive!(fg, 
                              lastPose::Symbol; 
                              refKey::Symbol=:simulated, 
                              solvable::Int=1,
                              graphinit::Bool=true,
                              landmarkRegex::Regex=r"l\d+",
                              srcNumber::Integer = match(r"\d+", string(lastPose)).match |> x->parse(Int,x),
                              atol::Real=1,
                              _doHoneycomb::Bool=false )
  #
  newFactor = RoME.Pose2Point2BearingRange(Normal(0,0.03), Normal(20,0.5))
  isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, lastPose, landmarkRegex, RoME.Point2, newFactor; atol=atol)

  # FIXME, oddball option that should be refactored (breaks Beehive usecase if true)
  if _doHoneycomb
    # force isAlready until fixed _checkVariableByReference parametric solution at -pi Optim issue
    global _honeycombRecipe
    isAlready = false
    genLabel = Symbol(:l, srcNumber)
    if haskey(_honeycombRecipe, genLabel) 
      isAlready = true
      genLabel = _honeycombRecipe[genLabel]
      simPPE = getPPE(fg, genLabel, refKey)
    end
  end

  # maybe add new variable
  if !isAlready
    @info "New variable with simPPE" genLabel round.(simPPE.suggested,digits=2)
    newVar = addVariable!(fg, genLabel, RoME.Point2, solvable=solvable, tags=[:LANDMARK;])
    addFactor!(fg, [lastPose; genLabel], newFactor, solvable=solvable, graphinit=graphinit)
    
    # also set :simulated PPE for similar future usage
    setPPE!(newVar, :simulated, typeof(simPPE), simPPE)   # TODO this API can be improved
  else
    @info "Adding simulated loop closure with perfect data association" lastPose genLabel
    addFactor!(fg, [lastPose; genLabel], newFactor, solvable=solvable, graphinit=graphinit)
  end

  #
  return genLabel
end


function _driveHex!(fgl::AbstractDFG,
                    posecount::Int;
                    poseCountTarget::Real=Inf,
                    graphinit::Bool=false,
                    gaugePrior::Symbol=:x0f1,
                    refKey::Symbol=:simulated,
                    addLandmarks::Bool=true,
                    landmarkSolvable::Int=1,
                    _doHoneycomb::Bool=false,
                    atol::Real=1,
                    postpose_cb::Function=(fg_,latestpose)->()  )
  #

  # Drive around in a hexagon
  for i in (posecount):(posecount+5)
    if poseCountTarget <= posecount
      break
    end
    psym = Symbol("x$i")
    # nsym = Symbol("x$(i+1)")
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    posecount += 1
    v_n = _addPoseCanonical!(fgl, psym, posecount, pp, refKey=refKey, graphinit=graphinit, postpose_cb=postpose_cb)
    nsym = getLabel(v_n)
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(fgl, nsym; refKey=refKey, solvable=landmarkSolvable, atol=atol, graphinit=false, _doHoneycomb=_doHoneycomb)
  end

  return posecount
end

function _offsetHexLeg( fgl::AbstractDFG,
                        posecount::Int; 
                        poseCountTarget::Real=Inf,
                        direction::Symbol=:right,
                        graphinit::Bool=false,
                        psym = Symbol("x$(posecount)"),
                        # nsym = Symbol("x$(posecount+1)"),
                        refKey::Symbol=:simulated,
                        guagePrior::Symbol=:x0f1,
                        addLandmarks::Bool=true,
                        landmarkSolvable::Int=1,
                        atol::Real=1,
                        _doHoneycomb::Bool=false,
                        postpose_cb::Function=(fg_,latestpose)->()   )
  #
  #skip out early if pose count target has been reached
  if poseCountTarget <= posecount
    return posecount
  end
  dirsign = if direction == :right
    -1
  elseif direction == :left
    +1
  else
    error("unknown direction symbol $direction")
  end
  pp = Pose2Pose2(MvNormal([10.0;0;dirsign*pi/3], diagm([0.1;0.1;0.1].^2)))
  
  posecount += 1
  v_n = _addPoseCanonical!(fgl, psym, posecount, pp, refKey=refKey, graphinit=graphinit, postpose_cb=postpose_cb)
  nsym = getLabel(v_n)

  # add a new landmark (if not yet present)
  !addLandmarks ? nothing : _addLandmarkBeehive!(fgl, nsym; refKey=refKey, solvable=landmarkSolvable, atol=atol, graphinit=false, _doHoneycomb=_doHoneycomb )

  return posecount
end


function generateGraph_Honeycomb!(poseCountTarget::Int=36;
                                  graphinit::Bool = false,
                                  dfg::AbstractDFG = LocalDFG{SolverParams}(solverParams=SolverParams(graphinit=graphinit)),  
                                  useMsgLikelihoods::Bool=getSolverParams(dfg).useMsgLikelihoods,
                                  direction::Symbol = :right,
                                  solvable::Int=1,
                                  refKey::Symbol=:simulated,
                                  addLandmarks::Bool=true,
                                  landmarkSolvable::Int=0,
                                  atol::Real=1,
                                  postpose_cb::Function=(fg_,latestpose)->()     )
  #
  global _honeycombRecipe

  # does anything exist in the graph yet
  posecount = if :x0 in ls(dfg)
    # what is the last pose
    lastPose = (ls(dfg, r"x\d+") |> sortDFG)[end]
    # get latest posecount number
    match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  else
    # initial zero pose
    generateGraph_ZeroPose(;dfg, varType=RoME.Pose2, postpose_cb, solverParams=SolverParams(;graphinit)) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue

    # # reference ppe on :x0
    # refVal = zeros(3)
    # ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
    # setPPE!(dfg[:x0], refKey, DFG.MeanMaxPPE, ppe)
    
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(dfg, :x0; refKey=refKey, solvable=landmarkSolvable, atol=atol, _doHoneycomb=true, graphinit=false)

    # staring posecount (i.e. :x0)
    0
  end

  # keep adding poses until the target number is reached
  while posecount < poseCountTarget
    posecount = _driveHex!(dfg, posecount; graphinit=graphinit, landmarkSolvable=landmarkSolvable, atol=atol, poseCountTarget=poseCountTarget, postpose_cb=postpose_cb, _doHoneycomb=true)
    # drive the offset legs
    lastPose = Symbol(:x, posecount)
    if haskey(_honeycombRecipe, lastPose)
      posecount = _offsetHexLeg(dfg, posecount; direction=_honeycombRecipe[lastPose], graphinit=graphinit, atol=atol, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget, postpose_cb=postpose_cb, _doHoneycomb=true)
    end
    posecount = _offsetHexLeg(dfg, posecount; direction=direction, graphinit=graphinit, atol=atol, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget, postpose_cb=postpose_cb, _doHoneycomb=true)
  end

  # NOTE solvable forced for everything at this time
  setSolvable!.(dfg, ls(dfg),  solvable)
  setSolvable!.(dfg, lsf(dfg), solvable)

  return dfg
end
