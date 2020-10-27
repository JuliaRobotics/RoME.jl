

using .Flux

import Base: convert
import IncrementalInference: getSample

export MixtureFluxPose2Pose2, PackedMixtureFluxPose2Pose2


struct MixtureFluxPose2Pose2{F <: FunctorInferenceType} <: AbstractRelativeRoots
  Zij::F
  # delta time between variables
  DT::Ref{Float64}
  specialSampler::Function
end

mutable struct PackedMixtureFluxPose2Pose2 <: PackedInferenceType
  packedZij::PackedMixture
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
                                outputDims::NTuple,
                                otherComponents::AbstractVector{<:SamplableBelief}=[MvNormal(zeros(3),diagm(ones(3)));],
                                diversity=[0.5;0.5];
                                DT::Real=0,
                                shuffle::Bool=true,
                                serializeHollow::Bool=false,
                                specialSampler::Function=sampleMixtureFluxPose2Pose2  )
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
  MixtureFluxPose2Pose2(mix, Base.RefValue{Float64}(DT), specialSampler)
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
  PackedMixtureFluxPose2Pose2(convert(Mixture,obj.packedZij),
                                      obj.DT,
                                      specFnc)
end




function calcVelocityInterPose2!( nfb::MixtureFluxPose2Pose2,
                                  iPts::AbstractMatrix{<:Real},
                                  jPts::AbstractMatrix{<:Real},
                                  idx::Int  )
  #
  joyVelData = nfb.Zij.components.fluxnn.data
  # DXY[1:2,i] .= TransformUtils.R(iPts[3,i])'*DXY[1:2,i]
  joyVelData[1,3:4] .= jPts[1:2,idx]
  joyVelData[1,3:4] .-= iPts[1:2,idx]
  joyVelData[1,3:4] ./= nfb.DT[] # convert to velocity
  joyVelData[1,3:4] .= TransformUtils.R(iPts[3,idx])'*joyVelData[1,3:4]
  # just set zero if something is wrong
  if isnan(joyVelData[1,3]) || isinf(abs(joyVelData[1,3])) || isnan(joyVelData[1,4]) || isinf(abs(joyVelData[1,4]))
    joyVelData[1,3:4] .= 0.0
  end
  for i in 2:size(joyVelData,1)
    joyVelData[i,3:4] .= joyVelData[1,3:4]
  end
end


function calcVelocityInterPose2!( nfb::MixtureFluxPose2Pose2,
                                  iPts::AbstractMatrix{<:Real},
                                  jPts::AbstractMatrix{<:Real}  )
  #
  @assert size(jPts,2) == size(iPts,2) "calcVelocityInterPose2! can currently only evaluate equal population size variables"

  # calculate an average velocity component
  if nfb.DT[] == 0
    # DXY = (@view jPts[1:2,:]) - (@view jPts[1:2,:])
    # rotate delta position from world to local iX frame
    for i in 1:size(iPts,2)
      calcVelocityInterPose2!(nfb, iPts, jPts, i)
    end
  else
    nfb.joyVelData[:,3:4] .= 0.0
  end
  nothing
end



function (nfb::MixtureFluxPose2Pose2)(
            res::AbstractArray{<:Real},
            fmd::FactorMetadata,
            idx::Int,
            meas::Tuple,
            Xi::AbstractArray{<:Real,2},
            Xj::AbstractArray{<:Real,2}  )
  #
  # if, use prediction sample
  if meas[2][idx] == 2
    # get live velocity estimate for each sample (nfb.joyVelData[:,3:4])
    calcVelocityInterPose2!(nfb, Xi, Xj, idx)
    # predict odom for this sample from a specific prediction model
    # NOTE, already integrated as part of Mixture
    # meas[1][1:2,idx] = nfb.allPredModels[meas[3][idx]](nfb.joyVelData)
  end

  # calculate the error for that measurement sample as Pose2Pose2
  nfb.Zij(res, fmd, idx, (meas[1],), Xi, Xj)
  nothing
end





## ============================================================================================
## Legacy support:  FluxModelsPose2Pose2
## ============================================================================================


export FluxModelsPose2Pose2


# Convenience function to help call the right constuctor
FluxModelsPose2Pose2( allModels::Vector{P},
                      jvd::D,
                      naiveModel::M,
                      naiveFrac::Real=0.5,
                      specialSampler::Function=sampleMixtureFluxPose2Pose2,
                      DT::Real=0.0,
                      shuffle::Bool=true ) where {P, M <: SamplableBelief, D <: AbstractMatrix} = MixtureFluxPose2Pose2(
                                        allModels,
                                        (25,4),
                                        jvd,
                                        (3,),
                                        [naiveModel;],
                                        [1-naiveFrac;naiveFrac],
                                        DT=DT,
                                        shuffle=shuffle,
                                        specialSampler=specialSampler )
#







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