

@info "RoME is adding Flux related functionality."

# Required packages
using .Flux

include("models/Pose2OdoNN_01.jl") # until a better way is found to deserialize

# make sure the lambdas get defined early
buildPose2OdoNN_01_FromElements();

# FluxModelsPose2Pose2

# the factor definitions
export FluxModelsPose2Pose2, sampleFluxModelsPose2Pose2, PackedFluxModelsPose2Pose2
# some utilities
export setShuffleAll!, setNaiveFracAll!



using Random, Statistics
using DistributedFactorGraphs, TransformUtils

import Base: convert
import IncrementalInference: getSample, calcZDim


struct FluxModelsPose2Pose2{P,D<:AbstractArray,M<:SamplableBelief} <: AbstractRelativeRoots
  allPredModels::Vector{P}
  joyVelData::D
  naiveModel::M
  naiveFrac::Ref{Float64}
  Zij::Pose2Pose2
  specialSampler::Function # special keyword field name used to invoke 'specialSampler' logic
  DT::Ref{Float64}
  shuffle::Ref{Bool}
end



function calcVelocityInterPose2!( nfb::FluxModelsPose2Pose2,
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


function calcVelocityInterPose2!( nfb::FluxModelsPose2Pose2,
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


function sampleFluxModelsPose2Pose2(nfb::FluxModelsPose2Pose2,
                                    N::Int,
                                    fmd::FactorMetadata,
                                    Xi::DFGVariable,
                                    Xj::DFGVariable )
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
  # can suppress shuffle for NN training purposes
  1 < numModels && nfb.shuffle[] ? shuffle!(allPreds) : nothing
  
  # cache the time difference estimate
  nfb.DT[] = (getTimestamp(Xj) - getTimestamp(Xi)).value * 1e-3
  
  return (smplsNaive, selectSource, allPreds)
end

# Convenience function to help call the right constuctor
FluxModelsPose2Pose2( allModels::Vector{P},
                      jvd::D,
                      naiveModel::M,
                      naiveFrac::Real=0.5,
                      ss::Function=sampleFluxModelsPose2Pose2,
                      DT::Real=0.0,
                      shuffle::Bool=true ) where {P, M <: SamplableBelief, D <: AbstractMatrix} = FluxModelsPose2Pose2{P,D,M}(
                                        allModels,
                                        jvd,
                                        naiveModel,
                                        naiveFrac,
                                        Pose2Pose2(MvNormal(zeros(3),diagm(ones(3)))), # this dummy distribution does not get used
                                        ss,
                                        DT,
                                        shuffle )
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

"""
    $SIGNATURES

Helper function to assemble the data required for each model.
"""
function assembleNNFactorData(dfg::AbstractDFG, fctsym::Symbol, idx::Int)
  #
  f1 = getFactor(dfg, fctsym)
  nfb = getFactorType(f1)
  varsyms = getVariableOrder(f1)
  Xi = getVariable(dfg, varsyms[1]) |> getBelief |> getPoints
  Xj = getVariable(dfg, varsyms[2]) |> getBelief |> getPoints

  # add the specific velocity information
  calcVelocityInterPose2!(nfb, Xi, Xj, idx)
  nfb.joyVelData
end


## packing converters

struct PackedFluxModelsPose2Pose2 <: IncrementalInference.PackedInferenceType
  joyVelData::Vector{Vector{Float64}}
  naiveModel::String
  naiveFrac::Float64
  modelWeightsAll::Vector{Dict{Symbol,Vector}}
end
    # W1::Matrix{Float64}
    # b1::Vector{Float64}
    # W2::Matrix{Float64}
    # b2::Vector{Float64}
    # W3::Matrix{Float64}
    # b3::Vector{Float64}
  # allPredModels::Vector{P}
  # joyVelData::D
  # naiveModel::M
  # naiveFrac::Ref{Float64}
  # Zij::Pose2Pose2
  # specialSampler::Function # special keyword field name used to invoke 'specialSampler' logic
  # DT::Ref{Float64}

function convert(::Type{FluxModelsPose2Pose2}, d::PackedFluxModelsPose2Pose2)
  MDS = []
  for mdw in d.modelWeightsAll
    # @show keys(mdw) |> collect
    mW1 = mdw[:W1]
    @cast W1c[row,col] := mW1[row][col]
    W1 = W1c |> Matrix{Float32}
    @assert W1[1,1] - mdw[:W1][1][1] |> abs < 1e-6 "Something about matrix vector order in FluxModelsPose2Pose2 converter broke."
    @assert W1[1,2] - mdw[:W1][1][2] |> abs < 1e-6 "Something about matrix vector order in FluxModelsPose2Pose2 converter broke."
    @assert W1[2,1] - mdw[:W1][2][1] |> abs < 1e-6 "Something about matrix vector order in FluxModelsPose2Pose2 converter broke."
    b1 = mdw[:b1] |> Vector{Float32}
    mW2 = (mdw[:W2])
    @cast W2c[row,col] := mW2[row][col]
    W2 = W2c |> Matrix{Float32}
    b2 = mdw[:b2] |> Vector{Float32}
    mW3 = (mdw[:W3])
    @cast W3c[row,col] := mW3[row][col]
    W3 = W3c |> Matrix{Float32}
    b3 = mdw[:b3] |> Vector{Float32}
    md = buildPose2OdoNN_01_FromElements(W1,b1,W2,b2,W3,b3)
    push!(MDS, md)
  end
  @cast joyVel[row,col] := d.joyVelData[row][col]
  joyVelM = joyVel |> Matrix{Float64}
  FluxModelsPose2Pose2(MDS, joyVelM, extractdistribution(d.naiveModel), d.naiveFrac )
end

# need to overwrite for specialSampler
# TODO rather consolidate with getDimension(::Pose2Pose2)
function calcZDim(usrfnc::FluxModelsPose2Pose2, Xi::Vector{<:DFGVariable})
  return 3
end

function convert(::Type{PackedFluxModelsPose2Pose2}, d::FluxModelsPose2Pose2)
  MDS = []
  for md in d.allPredModels
    mddict =Dict{Symbol, AbstractArray}()
    @cast W1c[row][col] := md[1].W1[row,col]
    @cast W2c[row][col] := md[5].W[row,col]
    @cast W3c[row][col] := md[6].W[row,col]
    mddict[:W1] = Vector{Float32}[vec(W1c[i]) for i in 1:length(W1c)]
    mddict[:b1] = md[1].b1
    mddict[:W2] = Vector{Float32}[vec(W2c[i]) for i in 1:length(W2c)]
    # mddict[:W2] = Vector{Vector{Float32}}(Vector{Float32}(W2c))
    mddict[:b2] = md[5].b
    mddict[:W3] = Vector{Float32}[vec(W3c[i]) for i in 1:length(W3c)]
    # mddict[:W3] = Vector{Vector{Float32}}(Vector{Float32}(W3c))
    mddict[:b3] = md[6].b
    push!(MDS, mddict)
  end
  joyVel = Vector{Float64}[d.joyVelData[i,:] for i in 1:size(d.joyVelData,1)]
  PackedFluxModelsPose2Pose2(joyVel, string(d.naiveModel), d.naiveFrac[], MDS)
end




## Utilities


function setShuffleAll!(dfg::AbstractDFG, shuf::Bool)
  fs = lsf(dfg, FluxModelsPose2Pose2)
  (x->getFactorType(dfg, x).shuffle[] = shuf).(fs)
  nothing
end

function setNaiveFracAll!(dfg::AbstractDFG, frac::Real)
  fs = lsf(dfg, FluxModelsPose2Pose2)
  (x->getFactorType(dfg, x).naiveFrac[] = frac).(fs)
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




#
