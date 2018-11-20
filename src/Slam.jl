

"""
$(TYPEDEF)
"""
mutable struct SLAMWrapper
  fg::IncrementalInference.FactorGraph
  tree::NothingUnion{IncrementalInference.BayesTree}
  lndmidx::Int
  lastposesym::Symbol
  lastpose::SE3
  lbl2usrid::Dict{Symbol,Int}
  usrid2lbl::Dict{Int, Symbol}
  SLAMWrapper() = new()
  SLAMWrapper(a,b,c) = new(a,b,c,:x1, SE3(0), Dict{Symbol,Int}(), Dict{Int,Symbol}()) # TODO to be deprecated
end


function addOdoFG!(slaml::SLAMWrapper, odo::Pose3Pose3;
                  N::Int=100, ready::Int=1,
                  saveusrid::Int=-1)
  #
  vprev = getVert(slaml.fg, slaml.lastposesym)
  # vprev, X, nextn = getLastPose(fgl)
  npnum = parse(Int,string(slaml.lastposesym)[2:end]) + 1
  nextn = Symbol("x$(npnum)")
  vnext = addNode!(slaml.fg, nextn, Pose2, N=N, ready=ready)
  # vnext = addNode!(slaml.fg, nextn, getVal(vprev) ⊕ odo, N=N, ready=ready, labels=["POSE"])
  slaml.lastposesym = nextn
  fact = addFactor!(slaml.fg, [vprev;vnext], odo)

  if saveusrid > -1
    slaml.lbl2usrid[nextn] = saveusrid
    slaml.usrid2lbl[saveusrid] = nextn
  end
  return vnext, fact
end


function addposeFG!(slaml::SLAMWrapper,
      constrs::Vector{IncrementalInference.FunctorInferenceType};
      N::Int=100,
      ready::Int=1,
      saveusrid::Int=-1   )
  #
  vprev = getVert(slaml.fg, slaml.lastposesym)

  npnum = parse(Int,string(slaml.lastposesym)[2:end]) + 1
  nextn = Symbol("x$(npnum)")
  # preinit
  vnext = nothing
  if !haskey(slaml.fg.IDs, nextn)
    vnext = addNode!(slaml.fg, nextn, Pose2, N=N, ready=ready)
    # vnext = addNode!(slaml.fg, nextn, getVal(vprev), N=N, ready=ready, labels=["POSE"])
  else
    vnext = getVert(slaml.fg, nextn) #, api=localapi # as optimization, assuming we already have latest vnest in slaml.fg
  end
  slaml.lastposesym = nextn

  addsubtype(fgl::FactorGraph, vprev, vnext, cc::IncrementalInference.FunctorPairwise) = addFactor!(fgl, [vprev;vnext], cc)
  addsubtype(fgl::FactorGraph, vprev, vnext, cc::IncrementalInference.FunctorSingleton) = addFactor!(fgl, [vnext], cc)

  facts = Graphs.ExVertex[]
  PP = BallTreeDensity[]
  for cns in constrs
    fa = addsubtype(slaml.fg, vprev, vnext, cns)
    push!(facts, fa)
  end

  # set node val from new constraints as init
  val = predictbelief(slaml.fg, vnext, facts, N=N)
  setVal!(vnext, val)
  IncrementalInference.dlapi.updatevertex!(slaml.fg, vnext)

  if saveusrid > -1
    slaml.lbl2usrid[nextn] = saveusrid
    slaml.usrid2lbl[saveusrid] = nextn
  end
  return vnext, facts
end





#
