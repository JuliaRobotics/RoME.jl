

import Base: convert
import IncrementalInference: getSample

struct MixtureFlux{N,F<:FunctorInferenceType,S,T}
  mixture::Mixture{N,F,S,T}
  # special keyword field name used to invoke 'specialSampler' logic
  specialSampler::Function 
end

const _IIFListTypes = Union{<:AbstractVector, <:Tuple, <:NTuple}

function MixtureFlux( F_::FunctorInferenceType, 
                      compList::_IIFListTypes, 
                      diversity::Union{<:AbstractVector, <: NTuple}) 
  #
  mix = Mixuture(F_, )
  MixtureFlux{N,F,S,T}()
end

MixtureFlux(::Type{F}, w...;kw...) where F <: FunctorInferenceType = MixtureFlux(F(LA.I), w...;kw...)


struct MixtureFluxPose2Pose2{N,S,T} <: AbstractRelativeRoots
  mixture::Mixture{N,F,S,T}
  # delta time between variables
  DT::Ref{Float64}
end

# const MixtureFluxPose2Pose2{S,T} = MixtureFluxPose2Pose2{2,S,T}


function sampleMixtureFluxPose2Pose2( nfb::MixtureFluxPose2Pose2,
                                      N::Int,
                                      fmd::FactorMetadata,
                                      legacy... )
  #
  # only calculate once after default value of 0
  if nfb.DT[] == 0
    # cache the time difference estimate
    nfb.DT[] = (getTimestamp(fmd.fullvariables[2]) - getTimestamp(fmd.fullvariables[1])).value * 1e-3
  end

  return getSample(nfb.Zij, N)
end



function MixtureFluxPose2Pose2(fluxmodels::AbstractVector,
                                data::AbstractArray;
                                naive::SamplableBelief=MvNormal(zeros(3),diagm(ones(3))),
                                diversity=[0.5;0.5],
                                DT::Real=0,
                                shuffle::Bool=true  )
  #
  mix = Mixture(Pose2Pose2, )
end







#