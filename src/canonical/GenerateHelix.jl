
export 
  generateCanonicalFG_Helix2D!,
  generateCanonicalFG_Helix2DSlew!,
  generateCanonicalFG_Helix2DSpiral!


"""
    $SIGNATURES

Generalized canonical graph generator function for helix patterns.

Notes
- assumes poses are labeled according to r"x\\d+"
- Gradient (i.e. angle) calculations are on the order of 1e-8.
- Use callback `spine_t(t)::Complex` to modify how the helix pattern is moved in x, y along the progression of `t`,
  - See related wrapper functions for convenient generators of helix patterns in 2D,
  - Real valued `xr_t(t)` and `yr_t(t)` can be modified (and will override) complex valued `spine_t` instead.
- use `postpose_cb = (fg_, lastestpose) -> ...` for additional user features after each new pose

Related

[`generateCanonicalFG_Helix2DSlew!`](@ref), [`generateCanonicalFG_Helix2DSpiral!`](@ref)
"""
function generateCanonicalFG_Helix2D!(numposes::Integer=40;
                                      posesperturn::Integer=20,
                                      graphinit::Bool=false,
                                      dfg::AbstractDFG = LightDFG{SolverParams}(solverParams=SolverParams(graphinit=graphinit)),
                                      useMsgLikelihoods::Bool=getSolverParams(dfg).useMsgLikelihoods,
                                      radius::Real=10,
                                      spine_t=(t)->0 + im*0,
                                      xr_t::Function=(t)->real(spine_t(t)),
                                      yr_t::Function=(t)->imag(spine_t(t)),
                                      poseRegex::Regex=r"x\d+",
                                      μ0=[0;0;pi/2],
                                      refKey::Symbol=:simulated,
                                      Qd::Matrix{<:Real}=diagm( [0.1;0.1;0.05].^2 ),
                                      postpose_cb::Function=(fg_,latestpose)->()   )
  #
  
  # add first pose if not already exists
  if !( :x0 in ls(dfg) )
    generateCanonicalFG_ZeroPose(dfg=dfg, μ0=μ0, graphinit=graphinit, postpose_cb=postpose_cb) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue
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
  tmp = calcHelix_T(0, turns, posesperturn, radius=radius, spine_t=spine_t, xr_t=xr_t, yr_t=yr_t)

  Tμ = SE2(μ0-[0;0;pi/2])

  bidx = 1
  eidx = 1
  for tn in 0:(ceil(Int, turns)-1)
    eidx += posesperturn
    eidx = minimum( [eidx, length(tmp[1])] )
    # tmp_ = _calcHelix2DApprox(N_ppt=posesperturn, radius=radius, runback=runback)
    tmp_ = hcat(tmp[2][bidx:eidx,:],tmp[3][bidx:eidx])'
    # adjust for turn progression in x
    # tmp_[1,:] .+= tn*(2radius*(1-runback))
    oldpose = Tμ*SE2(tmp_[:,1])
    
    # add each new pose (skippin the first element)
    for ps in 2:size(tmp_,2)
      # check exit condition
      numposes-1 <= posecount && break
      # add a new pose
      newpose = Tμ*TU.SE2(tmp_[:,ps])
      deltaodo = se2vee(oldpose \ newpose)
      factor = Pose2Pose2( MvNormal(deltaodo, Qd) )
      posecount += 1
      v_n = _addPoseCanonical!(dfg, lastpose, posecount, factor, poseRegex=poseRegex, refKey=refKey, overridePPE=se2vee(newpose), postpose_cb=postpose_cb)
      lastpose = getLabel(v_n)
      oldpose = newpose
    end

    bidx = eidx
  end
  
  # 
  return dfg
end


"""
    $SIGNATURES

Generate canonical slewed helix graph (like a flattened slinky).

Notes
- Use `slew_x` and `slew_y` to pull the "slinky" out in different directions at constant rate.
- See generalized helix generator for more details. 
- Defaults are choosen to slew along x and have multple trajectory intersects between consecutive loops of the helix.

Related

[`generateCanonicalFG_Helix2D!`](@ref), [`generateCanonicalFG_Helix2DSpiral!`](@ref)
"""
generateCanonicalFG_Helix2DSlew!( numposes::Integer=40;
                                  slew_x::Real=2/3,
                                  slew_y::Real=0,
                                  spine_t=(t)->slew_x*t + im*slew_y*t,
                                  kwargs...  ) = generateCanonicalFG_Helix2D!(numposes; spine_t=spine_t, kwargs...)
#

"""
    $SIGNATURES

Generate canonical helix graph that expands along a spiral pattern, analogous flower petals.

Notes
- This function wraps the complex `spine_t(t)` function to generate the spiral pattern.
  - `rate_a` and `rate_r` can be varied for different spiral behavior.
- See generalized helix generator for more details. 
- Defaults are choosen to slewto have multple trajectory intersects between consecutive loops of the helix and do a decent job of moving around coverage area with a relative balance of encircled area sizes.

Related 

[`generateCanonicalFG_Helix2D!`](@ref), [`generateCanonicalFG_Helix2DSlew!`](@ref)
"""
generateCanonicalFG_Helix2DSpiral!( numposes::Integer=100;
                                    rate_r=0.6,
                                    rate_a=6,
                                    spine_t=(t)->rate_r*(t^0.5)*cis(rate_a*(t^0.4)),
                                    kwargs...  ) = generateCanonicalFG_Helix2D!(numposes; spine_t=spine_t, kwargs...)
#


#