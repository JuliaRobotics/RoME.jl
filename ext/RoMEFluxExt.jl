module RoMEFluxExt

@info "RoME.jl is loading extension functionality related to Flux.jl."

using Flux

import Base: convert
import IncrementalInference: getSample

import RoME: MixtureFluxPose2Pose2, PackedMixtureFluxPose2Pose2, FluxModelsPose2Pose2


include("Pose2OdoNN_01.jl")


function getSample(cfo::CalcFactor{<:MixtureFluxPose2Pose2})
  #
  nfb = cfo.factor
  # fmd = cfo.metadata

  # only calculate once after default value of 0
  if nfb.DT[] == 0
    # cache the time difference estimate
    nfb.DT[] = (getTimestamp(cfo.fullvariables[2]) - getTimestamp(cfo.fullvariables[1])).value * 1e-3
  end

  cf_ = CalcFactor( cfo.factor.Zij, 0, cfo._legacyParams, cfo.cache, cfo.fullvariables, cfo.solvefor, cfo.manifold)

  smpl = getSample(cf_)

  return (smpl, cf_.factor.labels)

end


function MixtureFluxPose2Pose2( fluxmodels::AbstractVector,
                                inputDims::NTuple,
                                data::AbstractArray,
                                outputDims::NTuple,
                                otherComponents::AbstractVector{<:SamplableBelief}=[MvNormal(zeros(3),diagm(ones(3)));],
                                diversity=[0.5;0.5];
                                DT::Real=0,
                                shuffle::Bool=true,
                                serializeHollow::Bool=false )
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
  MixtureFluxPose2Pose2(mix, Base.RefValue{Float64}(DT))
end


function Base.convert(::Union{Type{<:AbstractPackedFactor},Type{<:PackedMixtureFluxPose2Pose2}},
                      obj::MixtureFluxPose2Pose2 )
  #
  PackedMixtureFluxPose2Pose2(convert(PackedMixture,obj.Z),
                              obj.DT[] )
end

function Base.convert(::Union{Type{<:AbstractFactor},Type{<:MixtureFluxPose2Pose2}},
                      obj::PackedMixtureFluxPose2Pose2 )
  #
  MixtureFluxPose2Pose2(convert(Mixture,obj.Z),
                                Ref(obj.DT) )
end




function calcVelocityInterPose2!( nfb::MixtureFluxPose2Pose2,
                                  iPts::AbstractVector{<:Real},
                                  jPts::AbstractVector{<:Real} )
  #
  joyVelData = nfb.Zij.components.fluxnn.data
  # DXY[1:2,i] .= TransformUtils.R(iPts[3,i])'*DXY[1:2,i]
  joyVelData[1,3:4] .= jPts[1:2]
  joyVelData[1,3:4] .-= iPts[1:2]
  joyVelData[1,3:4] ./= nfb.DT[] # convert to velocity
  joyVelData[1,3:4] .= TransformUtils.R(iPts[3])'*joyVelData[1,3:4]
  # just set zero if something is wrong
  if isnan(joyVelData[1,3]) || isinf(abs(joyVelData[1,3])) || isnan(joyVelData[1,4]) || isinf(abs(joyVelData[1,4]))
    joyVelData[1,3:4] .= 0.0
  end
  for i in 2:size(joyVelData,1)
    joyVelData[i,3:4] .= joyVelData[1,3:4]
  end
  return nothing
end


# function calcVelocityInterPose2!( nfb::MixtureFluxPose2Pose2,
#                                   iPts::AbstractMatrix{<:Real},
#                                   jPts::AbstractMatrix{<:Real}  )
#   #
#   @assert size(jPts,2) == size(iPts,2) "calcVelocityInterPose2! can currently only evaluate equal population size variables"

#   # calculate an average velocity component
#   if nfb.DT[] == 0
#     # DXY = (@view jPts[1:2,:]) - (@view jPts[1:2,:])
#     # rotate delta position from world to local iX frame
#     for i in 1:size(iPts,2)
#       calcVelocityInterPose2!(nfb, iPts, jPts, i)
#     end
#   else
#     nfb.joyVelData[:,3:4] .= 0.0
#   end
#   nothing
# end


function (cfo::CalcFactor{<:MixtureFluxPose2Pose2})(meas1, meas2, Xi, Xj)
  #
  # userdata = cfo.metadata
  # fmd = cfo.metadata  
  nfb = cfo.factor

  # if, use prediction sample
  if meas2 == 2
    # get live velocity estimate for each sample (nfb.joyVelData[:,3:4])
    calcVelocityInterPose2!(nfb, Xi, Xj)
    # predict odom for this sample from a specific prediction model
    # NOTE, already integrated as part of Mixture
    # meas[1][1:2,idx] = nfb.allPredModels[meas[3][idx]](nfb.joyVelData)
  end

  # calculate the error for that measurement sample as Pose2Pose2
  #TODO, this is constructor in the hot loop is way to expensive.
  cfZ = CalcFactor( cfo.factor.Z.mechanics, 0, cfo._legacyParams, cfo.cache, cfo.fullvariables, cfo.solvefor, cfo.manifold)

  return  cfZ(meas1, Xi, Xj)

end





## ============================================================================================
## Legacy support:  FluxModelsPose2Pose2
## ============================================================================================


# Convenience function to help call the right constuctor
FluxModelsPose2Pose2( 
  allModels::Vector{P},
  jvd::D,
  naiveModel::M,
  naiveFrac::Real=0.5,
  DT::Real=0.0,
  shuffle::Bool=true 
) where {P, M <: SamplableBelief, D <: AbstractMatrix} = MixtureFluxPose2Pose2(
    allModels,
    (25,4),
    jvd,
    (3,),
    [naiveModel;],
    [1-naiveFrac;naiveFrac],
    DT=DT,
    shuffle=shuffle 
  )
#







# struct MixtureFlux{N,F<:AbstractFactor,S,T}
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

# function MixtureFlux( F_::AbstractFactor, 
#                       compList::_IIFListTypes, 
#                       diversity::Union{<:AbstractVector, <:NTuple, <:DiscreteNonParametric}) 
#   #
#   mix = Mixture(F_, compList, diversity)
#   _populateMixture(mix{N,F,S,T}) where {N,F,S,T} = MixtureFlux{N,F,S,T}(mix, sampleMixtureFlux)
#   return _populateMixture(mix)
# end

# MixtureFlux(::Type{F}, w...;kw...) where F <: AbstractFactor = MixtureFlux(F(LA.I), w...;kw...)




#