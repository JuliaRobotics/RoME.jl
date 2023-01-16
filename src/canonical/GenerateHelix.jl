


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
- can be used to grow a graph with repeated calls, but keyword parameters are assumed identical between calls.

See also: [`generateGraph_Helix2DSlew!`](@ref), [`generateGraph_Helix2DSpiral!`](@ref), [`generateGraph_Beehive!`](@ref)
"""
function generateGraph_Helix2D!(numposes::Integer = 40;
                                posesperturn::Integer = 20,
                                graphinit = nothing,
                                useMsgLikelihoods = nothing,
                                solverParams::SolverParams = SolverParams(;graphinit=false),
                                dfg::AbstractDFG = LocalDFG{SolverParams}(;solverParams),
                                radius::Real = 10,
                                spine_t = (t)->0 + im*0,
                                xr_t::Function = (t)->real(spine_t(t)),
                                yr_t::Function = (t)->imag(spine_t(t)),
                                poseRegex::Regex = r"x\d+",
                                μ0 = [0;0;pi/2],
                                refKey::Symbol = :simulated,
                                Qd::AbstractMatrix{<:Real} = diagm( [0.1;0.1;0.05].^2 ),
                                postpose_cb::Function = (fg_,latestpose)->()   )
  #
  (graphinit isa Nothing) ? nothing : @error("generateGraph_Helix2D! keyword graphinit obsolete, use solverParams=SolverParams(graphinit=..) instead.")
  (useMsgLikelihoods isa Nothing) ? nothing : @error("generateGraph_Helix2D! keyword useMsgLikelihoods obsolete, use solverParams=SolverParams(useMsgLikelihoods=..) instead.")
  
  # add first pose if not already exists
  _initpose = Symbol(match(r"[A-Za-z]+", poseRegex.pattern).match, 0)
  if !exists( dfg, _initpose )
    generateGraph_ZeroPose(;dfg, μ0, solverParams, postpose_cb) # , μ0=[0;0;1e-5] # tried for fix NLsolve on wrap issue
    # getSolverParams(dfg).useMsgLikelihoods = useMsgLikelihoods    
    # reference ppe on :x0
    ppe = DFG.MeanMaxPPE(refKey, μ0, μ0, μ0)
    setPPE!(dfg[:x0], refKey, DFG.MeanMaxPPE, ppe)
  end
  
  # start from existsing poses
  _poses = ls(dfg, poseRegex) |> sortDFG
  # what is the last pose and posecount number
  lastpose = _poses[end]
  posecount = match(r"\d+", string(lastpose)).match |> x->parse(Int,x)    
  # init how many poses at the beginning of a new turn
  bidx = length(_poses) # 1
  
  # fractional number of turns to make in total (after all graph generation is done)
  turns = numposes/posesperturn
  # generate helix pattern algebraically
  tmp = calcHelix_T(0, turns, posesperturn, radius=radius, spine_t=spine_t, xr_t=xr_t, yr_t=yr_t)
  # TODO, dont always start from 0 -- i.e. chop first repeat elements from deterministic helix
  
  # select the starting point
  _μ0 = μ0
  # @show _μ0 = 1 == bidx ? μ0 : getPPE(dfg, lastpose, refKey).suggested
  Tμ = SE2(_μ0-[0;0;pi/2]) # TODO update to Manifolds.jl
  
  # current end pose count for number of turns
  eidx = 1
  for tn in 0:(ceil(Int, turns)-1)
    eidx += posesperturn
    # skip out early if extending a previous existing graph
    eidx = minimum( [eidx, length(tmp[1])] )
    # tmp_ = _calcHelix2DApprox(N_ppt=posesperturn, radius=radius, runback=runback)
    tmp_ = hcat(tmp[2][bidx:eidx,:],tmp[3][bidx:eidx])'
    # adjust for turn progression in x
    # tmp_[1,:] .+= tn*(2radius*(1-runback))
    oldpose = Tμ*SE2(tmp_[:,1])
    eidx < bidx ? continue : nothing
    
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

[`generateGraph_Helix2D!`](@ref), [`generateGraph_Helix2DSpiral!`](@ref)
"""
generateGraph_Helix2DSlew!( numposes::Integer=40;
                            slew_x::Real=2/3,
                            slew_y::Real=0,
                            spine_t=(t)->slew_x*t + im*slew_y*t,
                            kwargs...  ) = generateGraph_Helix2D!(numposes; spine_t=spine_t, kwargs...)
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

[`generateGraph_Helix2D!`](@ref), [`generateGraph_Helix2DSlew!`](@ref)
"""
generateGraph_Helix2DSpiral!( numposes::Integer=100;
                              rate_r=0.6,
                              rate_a=6,
                              spine_t=(t)->rate_r*(t^0.5)*cis(rate_a*(t^0.4)),
                              kwargs...  ) = generateGraph_Helix2D!(numposes; spine_t=spine_t, kwargs...)
#


#