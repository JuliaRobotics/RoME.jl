


typealias VoidUnion{T} Union{Void, T}

type SLAMWrapper
  fg::IncrementalInference.FactorGraph
  tree::VoidUnion{IncrementalInference.BayesTree}
  lndmidx::Int
  lastposesym::Symbol
  lastpose::SE3
  lbl2usrid::Dict{Symbol,Int64}
  usrid2lbl::Dict{Int64, Symbol}
  SLAMWrapper() = new()
  SLAMWrapper(a,b,c) = new(a,b,c,:x1, SE3(0), Dict{Symbol,Int64}(), Dict{Int64,Symbol}()) # TODO to be deprecated
end


function addOdoFG!(slaml::SLAMWrapper, odo::Pose3Pose3;
                  N::Int=100, ready::Int=1,
                  saveusrid::Int=-1)
  #
  vprev = getVert(slaml.fg, slaml.lastposesym)
  # vprev, X, nextn = getLastPose(fgl)
  npnum = parse(Int,string(slaml.lastposesym)[2:end]) + 1
  nextn = Symbol("x$(npnum)")
  vnext = addNode!(slaml.fg, nextn, getVal(vprev) ⊕ odo, N=N, ready=ready)
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
  vnext = addNode!(slaml.fg, nextn, getVal(vprev), N=N, ready=ready)
  slaml.lastposesym = nextn

  facts = Graphs.ExVertex[]
  PP = BallTreeDensity[]
  for cns = constrs
    fa = nothing
    if typeof(cns) == FunctorPairwise
      fa = addFactor!(slaml.fg, [vprev;vnext], cns)
    elseif typeof(cns) == FunctorSingleton
      fa = addFactor!(slaml.fg, [vnext], cns)
    else
      error("Unrecognized type $typeof(cns)")
    end
    # ppts = evalFactor2(slaml.fg, fa, vnext)
    # push!(PP, kde!(ppts))
    push!(facts, fa)
  end

  # set node val from new constraints as init
  val = predictbelief(slaml.fg, vnext, facts, N=N)

  if saveusrid > -1
    slaml.lbl2usrid[nextn] = saveusrid
    slaml.usrid2lbl[saveusrid] = nextn
  end
  return vnext, facts
end


#
