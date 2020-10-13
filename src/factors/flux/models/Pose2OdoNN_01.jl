

export buildPyNNModel_01_FromElements, buildPyNNModel_01_FromWeights, loadTfModelIntoFlux
export buildPose2OdoNN_01_FromElements, buildPose2OdoNN_01_FromWeights, loadPose2OdoNNModelIntoFlux


function buildPose2OdoNN_01_FromElements( W1::AbstractMatrix{<:Real}=zeros(4,8),
                                          b1::AbstractVector{<:Real}=zeros(8),
                                          W2::AbstractMatrix{<:Real}=zeros(8,48),
                                          b2::AbstractVector{<:Real}=zeros(8),
                                          W3::AbstractMatrix{<:Real}=zeros(2,8),
                                          b3::AbstractVector{<:Real}=zeros(2))
  #
  # W1 = randn(Float32, 4,8)
  # b1 = randn(Float32,8)
  modjl = Chain(
    x -> (x*W1)' .+ b1 .|> relu,
    x -> reshape(x', 25,8,1),
    x -> maxpool(x, PoolDims(x, 4)),
    x -> reshape(x[:,:,1]',:),
    Dense(48,8,relu),
    Dense(8,2)
  )

  modjl[5].W .= W2
  modjl[5].b .= b2

  modjl[6].W .= W3
  modjl[6].b .= b3

  return modjl
end

# As loaded from tensorflow get_weights
# Super specialized function
function buildPose2OdoNN_01_FromWeights(pywe)
  buildPose2OdoNN_01_FromElements(pywe[1], pywe[2][:], pywe[3]', pywe[4][:], pywe[5]', pywe[6][:])
end

# function buildPose2OdoNN_01_FromWeightsGPU(pywe)
#   buildPose2OdoNN_01_FromElements(collect(pywe[1]) |> gpu, collect(pywe[2][:]) |> gpu, collect(pywe[3]') |> gpu, collect(pywe[4][:]) |> gpu, collect(pywe[5]') |> gpu, collect(pywe[6][:]) |> gpu) |> gpu
# end

@deprecate buildPyNNModel_01_FromElements(args...) buildPose2OdoNN_01_FromElements(args...)
@deprecate buildPyNNModel_01_FromWeights(args...) buildPose2OdoNN_01_FromWeights(args...)


#
