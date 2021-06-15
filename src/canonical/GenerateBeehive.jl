# generate a randomized Beehive

export generateCanonicalFG_Beehive!


function generateCanonicalFG_Beehive!(poseCountTarget::Int=10;
                                      dfg::AbstractDFG = initfg(),
                                      direction::Symbol = :right,
                                      graphinit::Bool = false,
                                      solvable::Int=1,
                                      refKey::Symbol=:simulated,
                                      addLandmarks::Bool=true,
                                      landmarkSolvable::Int=0,
                                      useMsgLikelihoods::Bool=getSolverParams(dfg).useMsgLikelihoods,
                                      postpose_cb::Function=(fg_,latestpose)->()     )
  #

  # does anything exist in the graph yet
  posecount = if :x0 in ls(dfg)
    # what is the last pose
    lastPose = (ls(dfg, r"x\d+") |> sortDFG)[end]
    # get latest posecount number
    match(r"\d+", string(lastPose)).match |> x->parse(Int,x)
  else
    # pick a random staring direction
    yaw0 = ([0.0;-2pi/3;2pi/3])[rand(1:3)]
    μ0 = [0;0;yaw0]
    # initial zero pose
    generateCanonicalFG_ZeroPose(fg=dfg, varType=Pose2, μ0=μ0, graphinit=graphinit, postpose_cb=postpose_cb) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue
    
    # add a new landmark (if not yet present)
    !addLandmarks ? nothing : _addLandmarkBeehive!(dfg, :x0, refKey=refKey, solvable=landmarkSolvable, graphinit=false)

    # staring posecount (i.e. :x0)
    0
  end

  # keep adding poses until the target number is reached
  while posecount < poseCountTarget
    direction = rand(1:2) == 1 ? :left : :right
    posecount = _offsetHexLeg(dfg, posecount, direction=direction, graphinit=graphinit, 
                              landmarkSolvable=landmarkSolvable, poseCountTarget=poseCountTarget, 
                              postpose_cb=postpose_cb)
  end

  # NOTE solvable forced for everything at this time
  setSolvable!.(dfg, ls(dfg),  solvable)
  setSolvable!.(dfg, lsf(dfg), solvable)

  return dfg  
end


#