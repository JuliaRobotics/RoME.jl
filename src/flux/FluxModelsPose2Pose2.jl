# FluxModelsPose2Pose2

export FluxModelsPose2Pose2, sampleFluxModelsPose2Pose2, PackedFluxModelsPose2Pose2


using .Flux

using Random, Statistics
using DistributedFactorGraphs, TransformUtils

import Base: convert
import IncrementalInference: getSample


struct FluxModelsPose2Pose2{P,D<:AbstractArray,M<:SamplableBelief} <: FunctorPairwise
  allPredModels::Vector{P}
  joyVelData::D
  naiveModel::M
  naiveFrac::Ref{Float64}
  Zij::Pose2Pose2
  specialSampler::Function # special keyword field name used to invoke 'specialSampler' logic
  DT::Ref{Float64}
end



function calcVelocityInterPose2!(nfb::FluxModelsPose2Pose2,
                                 iPts::AbstractMatrix{<:Real},
                                 jPts::AbstractMatrix{<:Real},
                                 idx::Int  )
  #
  # DXY[1:2,i] .= TransformUtils.R(iPts[3,i])'*DXY[1:2,i]
  nfb.joyVelData[1,3:4] .= jPts[1:2,idx]
  nfb.joyVelData[1,3:4] .-= iPts[1:2,idx]
  nfb.joyVelData[1,3:4] ./= nfb.DT[] # convert to velocity
  nfb.joyVelData[1,3:4] .= TransformUtils.R(iPts[3,idx])'*nfb.joyVelData[1,3:4]
  # just set zero if something is wrong
  if isnan(nfb.joyVelData[1,3]) || isinf(abs(nfb.joyVelData[1,3])) || isnan(nfb.joyVelData[1,4]) || isinf(abs(nfb.joyVelData[1,4]))
    nfb.joyVelData[1,3:4] .= 0.0
  end
  for i in 2:size(nfb.joyVelData,1)
    nfb.joyVelData[i,3:4] .= nfb.joyVelData[1,3:4]
  end
end


function calcVelocityInterPose2!(nfb::FluxModelsPose2Pose2,
                                 iPts::AbstractMatrix{<:Real},
                                 jPts::AbstractMatrix{<:Real}  )
  #
  @assert size(jPts,2) == size(iPts,2) "sampleFluxModelsPose2Pose2 can currently only evaluate equal population size variables"

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


# function FMP2P2Other()
#   calcVelocityInterPose2!(nfb, getVal(Xi), getVal(Xj) )
#   # sample predictive fraction
#
#   # and predict
#       # A = [rand(4) for i in 1:25]
#   smplsPredAll = nfb.predictFnc(nfb.joyVelData)
#   # add rotation dimension if not done by the prediction function
#   if size(smplsPredAll,2) == 2
#     smplsPredAll = hcat(smplsPredAll, zeros(size(smplsPredAll,1)))
#   end
#
#
#   Nn = sum(selectSource .== 2)
#   # calculate desired number of predicted values
#   Np = N - Nn
#   # randomly select particles for prediction (with possible duplicates forwhen Np > size(iPts,2))
#   Npp = Np < numModels ? Np : numModels
#   Nnn = N - Npp
#   selPreds = @view selectPredModel[1:Npp] # TODO better in-place
#   smpls_p = smplsPredAll[selPreds,:] # sample per row??
#   smpls_p[:,3] .= smplsNaive[3,Nnn+1:N] # use naive delta theta at this time
#   # naive fraction
#   smpls_n = @view smplsNaive[:, 1:Nnn] # sample per column
#   # join and shuffle predicted odo values
#   shfSmpl = shuffle!(1:N |> collect)
#   smpls = hcat(smpls_n, smpls_p')[:,shfSmpl]
# end


function sampleFluxModelsPose2Pose2(nfb::FluxModelsPose2Pose2,
                                    N::Int,
                                    fmd::FactorMetadata,
                                    Xi::DFGVariable,
                                    Xj::DFGVariable)::Tuple
 #

  # get the naive samples
  # model samples (all for theta at this time)
  smplsNaive = rand(nfb.naiveModel, N)

  # calculate naive model and Predictive fraction of samples, respectively
  selectSource = rand(Categorical([nfb.naiveFrac[]; 1-nfb.naiveFrac[]]), N)

  # number of predictors to choose from, and choose random subset
  numModels = length(nfb.allPredModels)
  allPreds = 1:numModels |> collect # 1:Npreds |> collect
  # TODO -- compensate when there arent enough prediction models
  if numModels < N
    repeat(allPreds, ceil(Int, (N-numModels)/N) + 1)
    allPreds = allPreds[1:N]
  end
  # samples for the order in which to use models, dont shuffle if N models
  1 < numModels && N != numModels ? shuffle!(allPreds) : nothing

  # cache the time difference estimate
  nfb.DT[] = (getTimestamp(Xj) - getTimestamp(Xi)).value * 1e-3

  return (smplsNaive, selectSource, allPreds)
end

# Convenience function to help call the right constuctor
FluxModelsPose2Pose2(allModels::Vector{P},
                     jvd::D,
                     naiveModel::M,
                     naiveFrac::Real=0.5,
                     ss::Function=sampleFluxModelsPose2Pose2,
                     DT::Real=0.0  ) where {P, M <: SamplableBelief, D <: AbstractMatrix} = FluxModelsPose2Pose2{P,D,M}(
                                        allModels,
                                        jvd,
                                        naiveModel,
                                        naiveFrac,
                                        Pose2Pose2(MvNormal(zeros(3),diagm(ones(3)))), # this dummy distribution does not get used
                                        ss,
                                        DT  )
#


function (nfb::FluxModelsPose2Pose2)(
            res::AbstractArray{<:Real},
            userdata::FactorMetadata,
            idx::Int,
            meas::Tuple{AbstractArray{<:Real},AbstractArray{<:Int},AbstractArray{<:Int}},
            Xi::AbstractArray{<:Real,2},
            Xj::AbstractArray{<:Real,2}  )
  #
  # if, use prediction sample
  if meas[2][idx] == 2
    # get live velocity estimate for each sample (nfb.joyVelData[:,3:4])
    calcVelocityInterPose2!(nfb, Xi, Xj, idx)
    # predict odom for this sample from a specific prediction model
    meas[1][1:2,idx] = nfb.allPredModels[meas[3][idx]](nfb.joyVelData)
  end

  # calculate the error for that measurement sample as Pose2Pose2
  nfb.Zij(res, userdata, idx, (meas[1],), Xi, Xj)
  nothing
end



## packing converters

struct PackedFluxModelsPose2Pose2 <: IncrementalInference.PackedInferenceType
  joyVelData::Matrix{Float64}
  naiveModel::String
  naiveFrac::Float64
end


function convert(::Type{FluxModelsPose2Pose2}, d::PackedFluxModelsPose2Pose2)
  FluxModelsPose2Pose2(PyTFOdoPredictorPoint2,d.joyVelData,extractdistribution(d.naiveModel),d.naiveFrac)
end

function convert(::Type{PackedFluxModelsPose2Pose2}, d::FluxModelsPose2Pose2)
  PackedFluxModelsPose2Pose2(d.joyVelData, string(d.naiveModel), d.naiveFrac)
end



#
