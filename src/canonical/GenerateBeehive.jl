

global _beehiveRecipe = Dict{Symbol,Symbol}(
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
                              graphinit::Bool=true  )
  #
  newFactor = RoME.Pose2Point2BearingRange(Normal(0,0.03), Normal(20,0.5))
  isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, lastPose, r"l\\d+", Point2, newFactor)

  # force isAlready until fixed _checkVariableByReference parametric solution at -pi Optim issue
  global _beehiveRecipe
  isAlready = false
  srcNumber = match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  genLabel = Symbol(:l, srcNumber)
  if haskey(_beehiveRecipe, genLabel) 
    isAlready = true
    genLabel = _beehiveRecipe[genLabel]
    simPPE = getPPE(fg, genLabel, refKey)
  end

  # maybe add new variable
  if !isAlready
    @info "New variable with simPPE" genLabel round.(simPPE.suggested,digits=2)
    newVar = addVariable!(fg, genLabel, Point2, solvable=solvable)
    addFactor!(fg, [lastPose; genLabel], newFactor, solvable=solvable, graphinit=graphinit)
    
    # also set :simulated PPE for similar future usage
    # newPPE = DFG.MeanMaxPPE(:simulated, simPPE, simPPE, simPPE)
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
    v_n = _addPose2Canonical!(fgl, psym, posecount, pp, refKey=refKey, graphinit=graphinit, postpose_cb=postpose_cb)
    nsym = getLabel(v_n)
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(fgl, nsym, refKey=refKey, solvable=landmarkSolvable, graphinit=false)
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
  v_n = _addPose2Canonical!(fgl, psym, posecount, pp, refKey=refKey, graphinit=graphinit, postpose_cb=postpose_cb)
  nsym = getLabel(v_n)

  # add a new landmark (if not yet present)
  !addLandmarks ? nothing : _addLandmarkBeehive!(fgl, nsym, refKey=refKey, solvable=landmarkSolvable, graphinit=false)

  return posecount
end


function generateCanonicalFG_Beehive!(poseCountTarget::Int=36;
                                      fg::AbstractDFG = initfg(),
                                      direction::Symbol = :right,
                                      graphinit::Bool = false,
                                      refKey::Symbol=:simulated,
                                      addLandmarks::Bool=true,
                                      landmarkSolvable::Int=0,
                                      useMsgLikelihoods::Bool=false     )
  #
  global _beehiveRecipe

  # does anything exist in the graph yet
  posecount = if :x0 in ls(fg)
    # what is the last pose
    lastPose = (ls(fg, r"x\d+") |> sortDFG)[end]
    # get latest posecount number
    match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  else
    # initial zero pose
    generateCanonicalFG_ZeroPose2(fg=fg, graphinit=graphinit) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue
    getSolverParams(fg).useMsgLikelihoods = useMsgLikelihoods    

    # reference ppe on :x0
    refVal = zeros(3)
    ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
    setPPE!(fg[:x0], refKey, DFG.MeanMaxPPE, ppe)
    
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(fg, :x0, refKey=refKey, solvable=landmarkSolvable, graphinit=false)

    # staring posecount (i.e. :x0)
    0
  end

  # keep adding poses until the target number is reached
  while posecount < poseCountTarget
    posecount = _driveHex!(fg, posecount, graphinit=graphinit, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget)
    # drive the offset legs
    lastPose = Symbol(:x, posecount)
    if haskey(_beehiveRecipe, lastPose)
      posecount = _offsetHexLeg(fg, posecount, direction=_beehiveRecipe[lastPose], graphinit=graphinit, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget)    
    end
    posecount = _offsetHexLeg(fg, posecount, direction=direction, graphinit=graphinit, landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget)
  end

  # FIXME force everything solvable for now
  setSolvable!.(fg, ls(fg), 1)
  setSolvable!.(fg, lsf(fg), 1)

  return fg
end
