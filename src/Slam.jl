


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
