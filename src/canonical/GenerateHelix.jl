
export generateCanonicalFG_Helix2D!
export generateCanonicalFG_Helix2DSlew!, generateCanonicalFG_Helix2DSpiral!



"""
    $SIGNATURES

Generate generalized helix parameterized by a curve along "t-axis" (i.e. z-axis, assuming z(t)=t).  

Notes
- Returns vectors for (`t`, `x,y`, and `yaw` angle).
- Offset to start at origin and facing direction along +y-axis.
- Use callbacks `xr_t(t)` and `yr_t(t)` to skew the helix with any desired curve, examples include
  - `xr_t = (t) -> (1/3)t` to generate helix pattern along x-axis,
  - or make spiral along t using xr_t, yr_t to generate a rose pattern on xy,
  - use `spine_t(t)=xr_t(t) + im*yr_t(t)` as shortcut for more complicated patterns,
  - note `xr_t` and `yr_t` are scaled by a factor `radius`, unscale the input by division if desired.
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
                      spine_t=(t)->0 + im*0,
                      xr_t::Function=(t)->real(spine_t(t)),
                      yr_t::Function=(t)->imag(spine_t(t))  )
  #
  # calc the position
  f(t, x=xr_t(t), y=yr_t(t)) = radius*( cis(pi - 2pi*t) + 1 + x + im*y)
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
                                      # runback::Real=2/3,
                                      spine_t=(t)->0 + im*0,
                                      xr_t::Function=(t)->real(spine_t(t)),
                                      yr_t::Function=(t)->imag(spine_t(t)),
                                      graphinit::Bool=false,
                                      poseRegex::Regex=r"x\d+",
                                      useMsgLikelihoods::Bool=true,
                                      refKey::Symbol=:simulated,
                                      Qd::Matrix{<:Real}=diagm( [0.1;0.1;0.05].^2 ),
                                      postpose_cb::Function=(fg_,latestpose)->()   )
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
  lastpose = (ls(dfg, poseRegex) |> sortDFG)[end]
  # get latest posecount number
  posecount = match(r"\d+", string(lastpose)).match |> x->parse(Int,x)    
  
  turns = numposes/posesperturn
  # TODO dont always start from 0
  tmp = _calcHelix_T(0, turns, posesperturn, radius=radius, spine_t=spine_t, xr_t=xr_t, yr_t=yr_t)

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
      v_n = _addPose2Canonical!(dfg, lastpose, posecount, factor, poseRegex=poseRegex, refKey=refKey, overridePPE=tmp_[:,ps], postpose_cb=postpose_cb)
      lastPose = getLabel(v_n)
      oldpose = newpose
    end

    bidx = eidx
  end
  
  # 
  return dfg
end



generateCanonicalFG_Helix2DSlew!( numposes::Integer=40;
                                  slew_x::Real=2/3,
                                  slew_y::Real=0,
                                  spine_t=(t)->slew_x*t + im*slew_y*t,
                                  kwargs...  ) = generateCanonicalFG_Helix2D!(numposes; spine_t=spine_t, kwargs...)
#


generateCanonicalFG_Helix2DSpiral!( numposes::Integer=100;
                                    rate_r=0.3,
                                    rate_a=4,
                                    spine_t=(t)->rate_r*(t^0.5)*cis(rate_a*(t^0.4)),
                                    kwargs...  ) = generateCanonicalFG_Helix2D!(numposes; spine_t=spine_t, kwargs...)
#


#