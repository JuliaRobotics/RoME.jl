
export Pose2Point2, PackedPose2Point2

#-------------------------------------------------------------------------------
#

"""
    $TYPEDEF

Bearing and Range constraint from a Pose2 to Point2 variable.
"""
mutable struct Pose2Point2{T <: IIF.SamplableBelief} <: IncrementalInference.FunctorPairwiseMinimize
    Zij::T
    Pose2Point2{T}() where {T} = new{T}()
    Pose2Point2{T}(x1::T) where {T <: IIF.SamplableBelief} = new{T}(x1)
end
Pose2Point2(x1::T) where {T <: IIF.SamplableBelief} = Pose2Point2{T}(x1)
function getSample(pp2br::Pose2Point2, N::Int=1)
  smpls = zeros(2, N)
  smpls[1:2,:] .= rand(pp2br.Zij, N)

  return (smpls,)
end
# define the conditional probability constraint
function (pp2br::Pose2Point2)(res::Array{Float64},
                              userdata::FactorMetadata,
                              idx::Int,
                              meas::Tuple{Array{Float64,2}},
                              xi::Array{Float64,2},
                              lm::Array{Float64,2} )
  #

  wLi = SE2(xi[:,idx])*SE2([meas[1][:,idx];0.0])
  res[1:2] .= lm[:,idx] .- se2vee(wLi)[1:2]

  res .*= res
  res[1] += res[2]
  res[2] = 0.0
  return res[1]
end


#TODO wrapper
# function (s::Pose2Point2{<:Normal})(xi::AbstractVector{T}, lm::AbstractVector{T}; kwargs...) where T <: Real
#   meas = [mean(s.bearing), mean(s.range)]
#   iΣ = [var(s.bearing)         0.0;
#                    0.0 var(s.range)]
#   #
#
#   # world frame
#   θ = meas[1] + xi[3]
#   mx = meas[2]*cos(θ)
#   my = meas[2]*sin(θ)
#
#   ex = lm[1] - (mx + xi[1])
#   ey = lm[2] - (my + xi[2])
#   er = sqrt(ex^2 + ey^2)
#
#   eθ = atan((my + xi[2]), (mx + xi[1])) - atan(lm[2], lm[1])
#
#   res = [eθ, er]
#
#   return res' * iΣ * res
# end


## Serialization support

mutable struct PackedPose2Point2 <: IncrementalInference.PackedInferenceType
    Zij::String
    PackedPose2Point2() = new()
    PackedPose2Point2(s1::AS) where {AS <: AbstractString} = new(string(s1))
end

function convert(::Type{PackedPose2Point2}, obj::Pose2Point2{T}) where {T <: IIF.SamplableBelief}
  return PackedPose2Point2(string(obj.Zij))
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert(::Type{Pose2Point2}, packed::PackedPose2Point2)
  Pose2Point2(extractdistribution(packed.Zij))
end
