

import Base: convert
import IncrementalInference: getSample

export MixtureFluxPose2Pose2. PackedMixtureFluxPose2Pose2


struct MixtureFluxPose2Pose2{F <: FunctorInferenceType} <: AbstractRelativeRoots
  Zij::F
  # delta time between variables
  DT::Ref{Float64}
  specialSampler::Function
end

mutable struct PackedMixtureFluxPose2Pose2
  packedMixture::PackedMixture
  DT::Float64
  specialSampler::String
end

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



function MixtureFluxPose2Pose2( fluxmodels::AbstractVector,
                                inputDims::NTuple,
                                data::AbstractArray,
                                outputDims::NTuple;
                                otherComponents::AbstractVector{<:SamplableBelief}=[MvNormal(zeros(3),diagm(ones(3)));],
                                diversity=[0.5;0.5],
                                DT::Real=0,
                                shuffle::Bool=true,
                                serializeHollow::Bool=false  )
  #
  mix = MixtureFluxModels(Pose2Pose2, 
                          fluxmodels, 
                          inputDims, 
                          data, 
                          outputDims, 
                          otherComponents, 
                          diversity, 
                          shuffle=shuffle, 
                          serializeHollow=serializeHollow)
  #
  MixtureFluxPose2Pose2(mix, DT, sampleMixtureFluxPose2Pose2)
end




function Base.convert(::Union{Type{<:PackedInferenceType},Type{<:PackedMixtureFluxPose2Pose2}},
                      obj::MixtureFluxPose2Pose2 )
  #
  toFnc = typeof(obj.specialSampler)
  fncName = string(toFnc.name.module)*"."*string(toFnc.name.name)
  PackedMixtureFluxPose2Pose2(convert(PackedMixture,obj.Zij),
                              obj.DT[],
                              fncName)
end

function Base.convert(::Union{Type{<:FunctorInferenceType},Type{<:MixtureFluxPose2Pose2}},
                      obj::PackedMixtureFluxPose2Pose2 )
  #
  sFnc = split(obj.specialSampler, '.') .|> Symbol
  specFnc = getfield(getfield(Main, sFnc[1]), sFnc[2])
  PackedMixtureFluxPose2Pose2(convert(Mixture,obj.packedMixture),
                                      obj.DT,
                                      specFnc)
end




# struct MixtureFlux{N,F<:FunctorInferenceType,S,T}
#   mixture::Mixture{N,F,S,T}
#   # special keyword field name used to invoke 'specialSampler' logic
#   specialSampler::Function 
# end

# # to be overrided by user via more specialized dispatch
# # FIXME as part of 
# function sampleMixtureFlux( nfb::MixtureFluxPose2Pose2,
#                             N::Int,
#                             fmd::FactorMetadata,
#                             legacy... )
#   #
#   getSample(nfb.mixture, N)
# end

# # const _IIFListTypes = Union{<:AbstractVector, <:Tuple, <:NTuple}

# function MixtureFlux( F_::FunctorInferenceType, 
#                       compList::_IIFListTypes, 
#                       diversity::Union{<:AbstractVector, <:NTuple, <:DiscreteNonParametric}) 
#   #
#   mix = Mixture(F_, compList, diversity)
#   _populateMixture(mix{N,F,S,T}) where {N,F,S,T} = MixtureFlux{N,F,S,T}(mix, sampleMixtureFlux)
#   return _populateMixture(mix)
# end

# MixtureFlux(::Type{F}, w...;kw...) where F <: FunctorInferenceType = MixtureFlux(F(LA.I), w...;kw...)




#