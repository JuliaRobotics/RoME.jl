# generate 2D boxy canonical graph pattern



_eps_w_int(::Type{T}) where {T <: Real} = eps(T)
_eps_w_int(::Type{T}) where {T <: Integer} = 1e-17


function _calcOdoBox( t4::T;
                      slew_x::Real=2/3,
                      length_x::Real = 15,
                      length_y::Real = length_x ) where {T <: Real}
  #

  # select direction
  t_ = t4 % 4
  if isapprox(t_, 0, atol=_eps_w_int(T))
    # drive along x
    return [length_x;0.0], :POSITIVE_X   
  elseif isapprox(t_, 1, atol=_eps_w_int(T))
    # drive along y
    return [0.0;length_y], :POSITIVE_Y
  elseif isapprox(t_, 2, atol=_eps_w_int(T))
    # drive along -x
    return [-slew_x*length_x;0.0], :NEGATIVE_X
  elseif isapprox(t_, 3, atol=_eps_w_int(T))
    # drive along -y    
    return [0.0;-length_y], :NEGATIVE_Y
  else
    error("non integer values of ")
  end
end


function driveLeg!( fg, 
                    lastPose,      
                    odo::AbstractVector{<:Real},
                    direction::Symbol;
                    graphinit::Bool=true,
                    Qd::AbstractArray{<:Real}=[1.0; 1.0],
                    factor::AbstractRelative=Point2Point2(MvNormal(odo, Qd)),
                    overridePPE=nothing,
                    postpose_cb::Function=(fg_,latestpose)->() )
  #
  # fg - factor graph object
  # lastpose
  # start - 
  # v - displacement vector (e.g. [length_x; 0.0] for direction= :NORTH) ASSUMES ONE ELEMENT ALWAYS ZERO
  # direction - direction symbol (:NORTH,:SOUTH,:EAST,:WEST)

  # add end pose of the leg
  newPose = incrSuffix(lastPose)
  
  v_n = RoME._addPoseCanonical!(fg, lastPose, -1, factor, genLabel=newPose, srcType=RoME.Point2,
                                graphinit=false, variableTags=[:POSE;direction], 
                                factorTags=[:ODOMETRY; direction],
                                overridePPE=overridePPE,
                                postpose_cb=postpose_cb  )
  #

  #
  # addFactor!(fg, [newPose;], PartialPrior(Normal(trueY, 0.1),(2,)), graphinit=graphinit)

  # return new pose name
  return v_n
end




# drive clockwise, x is North, y is East (NED convention).
# boxes start bottom left, spine of boxy helix is on x-axis
function driveOneBox!(fg;
                      lastPose=sortDFG(ls(fg, tags=[:POSE]))[end],
                      slew_x = 2/3,
                      length_x = 15,
                      length_y = length_x,
                      postpose_cb::Function=(fg_,latestpose)->()  )
  #
  # (g,lp) -> callbackSequenceScalar_ex4(g, lp, start_x, v )
  
  for leg in 0:3
    odo, direction = _calcOdoBox(leg, slew_x=slew_x, length_x=length_x, length_y=length_y)
    lastPose = driveLeg!(fg, lastPose, odo, direction, 
                          postpose_cb=postpose_cb  ) |> getLabel
    #
  end

  #
  nothing
end



"""
    $SIGNATURES

Canonical graph generator for box shapes along x and y axes, with default behaviour slewing each new box along x-axes with 2/3 overlap on edges.

Notes
- Function is still experimental
- TODO currently `numposes` works in multiples of 4, for each corner of a new box -- i.e. 8 --> 2 boxes. 
- Use `postpose_cb(fg, lastpose)` for additional behaviour after each pose variable has been added to the graph.

Related

[`generateGraph_ZeroPose`](@ref), [`generateGraph_Helix2DSlew!`](@ref), [`generateGraph_Hexagonal`](@ref)
"""
function generateGraph_Boxes2D!(numposes::Integer=16;
                                graphinit::Bool=false,
                                dfg::AbstractDFG = LocalDFG{SolverParams}(solverParams=SolverParams(graphinit=graphinit)),
                                useMsgLikelihoods::Bool=getSolverParams(dfg).useMsgLikelihoods,
                                length_x::Real=15,
                                length_y::Real=length_x,
                                slew_x::Real=2/3,
                                poseRegex::Regex=r"x\d+",
                                refKey::Symbol=:simulated,
                                Qd::Matrix{<:Real}=diagm( [0.1;0.1;0.05].^2 ),
                                postpose_cb::Function=(fg_,latestpose)->()   )
  #

  # actually start adding nodes 
  generateGraph_ZeroPose(dfg=dfg, varType=RoME.Point2, variableTags=[:POSE;], postpose_cb=postpose_cb)

  numboxes = ceil(Int, numposes/4)
  
  for tr in 0:(numboxes-1)
    # drive boxes
    driveOneBox!(dfg, slew_x=slew_x, length_x=length_x, postpose_cb=postpose_cb )
  end

  return dfg
end



#
