
export generateCanonicalFG_Helix2D!


"""
    $SIGNATURES

Generate generalized helix parameterized by a curve along "t-axis" (i.e. z-axis, assuming z(t)=t).  

Notes
- Returns vectors for (`t`, `x,y`, and `yaw` angle).
- Offset to start at origin and facing direction along +y-axis.
- Use callbacks `x_t(t)` and `y_t(t)` to skew the helix with any desired curve, examples include
  - `x_t = (t) -> (1/3)t` to generate helix pattern along x-axis,
  - or make spiral along t using x_t, y_t to generate a rose pattern on xy.
- Use the function twice for simulated and noisy trajectories (i.e. easier Gauss-Markov processes)
- Gradient (i.e. angle) calculations are on the order of 1e-8.

Related

[`generateCanonicalFG_Helix2D!`](@ref)
"""
function _calcHelix_T(start::Real=0,
                      stop::Real=1,
                      pointsperturn=20;
                      T::AbstractVector{<:Real}=(start:(stop*pointsperturn))./pointsperturn,
                      radius::Real = 0.5,
                      x_t::Function=(t)->0,
                      y_t::Function=(t)->0  )
  #
  # calc the position
  f(t, x=x_t(t), y=y_t(t)) = radius*( exp(im*(pi - 2pi*t)) + 1) + x + im*y
  vals = f.(T)

  # calc the gradient
  g(t, h=1e-8) = (f(t+h)-f(t))/h
  grad = g.(T)

  return T, hcat(real.(vals), imag.(vals)), angle.(grad)
end




## ================================================================================================
## GENERATE CANONICAL GRAPH
## ================================================================================================



# assume poses are labeled according to r"x\d+"
#- Gradient (i.e. angle) calculations are on the order of 1e-8.
function generateCanonicalFG_Helix2D!(numposes::Integer=40;
                                      posesperturn::Integer=20,
                                      dfg::AbstractDFG=initfg(),
                                      radius::Real=10,
                                      runback::Real=1/3,
                                      graphinit::Bool=false,
                                      poseRegex::Regex=r"x\d+",
                                      useMsgLikelihoods::Bool=true,
                                      refKey::Symbol=:simulated,
                                      Qd::Matrix{<:Real}=diagm( [0.1;0.1;0.05].^2 )   )
  #
  
  # add first pose if not already exists
  if !( :x0 in ls(dfg) )
    μ0=[0;0;pi/2]
    generateCanonicalFG_ZeroPose2(fg=dfg, μ0=μ0, graphinit=graphinit) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue
    getSolverParams(dfg).useMsgLikelihoods = useMsgLikelihoods    
    # reference ppe on :x0
    ppe = DFG.MeanMaxPPE(refKey, μ0, μ0, μ0)
    setPPE!(dfg[:x0], refKey, DFG.MeanMaxPPE, ppe)
  end
  
  # what is the last pose
  @show ls(dfg)
  lastPose = (ls(dfg, poseRegex) |> sortDFG)[end]
  # get latest posecount number
  posecount = match(r"\d+", string(lastPose)).match |> x->parse(Int,x)    
  
  turns = numposes/posesperturn
  # TODO dont always start from 0
  tmp = _calcHelix_T(0, turns, posesperturn, radius=radius, x_t=t->radius*runback*t)

  bidx = 1
  eidx = 1
  for tn in 0:(ceil(Int, turns)-1)
    eidx += posesperturn
    eidx = minimum( [eidx, length(tmp[1])] )
    # tmp_ = _calcHelix2DApprox(N_ppt=posesperturn, radius=radius, runback=runback)
    tmp_ = hcat(tmp[2][bidx:eidx,:],tmp[3][bidx:eidx])'
    # adjust for turn progression in x
    # tmp_[1,:] .+= tn*(2radius*(1-runback))
    oldpose = SE2(tmp_[:,1])
    
    # add each new pose (skippin the first element)
    for ps in 2:size(tmp_,2)
      # check exit condition
      numposes-1 <= posecount && break
      # add a new pose
      newpose = TU.SE2(tmp_[:,ps])
      deltaodo = se2vee(oldpose \ newpose)
      factor = Pose2Pose2( MvNormal(deltaodo, Qd) )
      posecount += 1
      v_n = _addPose2Canonical!(dfg, lastPose, posecount, factor, poseRegex=poseRegex, refKey=refKey, overridePPE=tmp_[:,ps])
      lastPose = getLabel(v_n)
      oldpose = newpose
    end

    bidx = eidx
  end
  
  # 
  return dfg
end




#