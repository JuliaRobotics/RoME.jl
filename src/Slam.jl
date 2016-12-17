


typealias VoidUnion{T} Union{Void, T}

type SLAMWrapper
  fg::IncrementalInference.FactorGraph
  tree::VoidUnion{IncrementalInference.BayesTree}
  lndmidx::Int
  lastposesym::Symbol
  SLAMWrapper() = new()
  SLAMWrapper(a,b,c) = new(a,b,c,:x0) # TODO to be deprecated
end


function addOdoFG!(slaml::SLAMWrapper, odo::Pose3Pose3;
                  N::Int=100, ready::Int=1)
    vprev, X, nextn = getLastPose(fgl)
    vnext = addNode!(fgl, nextn, X⊕odo, N=N, ready=ready)
    fact = addFactor!(fgl, [vprev;vnext], odo)
    return vnext, fact
end
