
struct P2P2BearingReuse
  measvec::Vector{Float64}
  predvec::Vector{Float64}
  resid::Vector{Float64}
  P2P2BearingReuse() = new(zeros(2),zeros(2),zeros(2))
end

"""
    $TYPEDEF

Single dimension bearing constraint from Pose2 to Point2 variable.
"""
struct Pose2Point2Bearing{B <: IIF.SamplableBelief} <: IncrementalInference.FunctorPairwiseMinimize
    bearing::B
    reuse::Vector{P2P2BearingReuse}
    Pose2Point2Bearing{B}() where B = new{B}()
    Pose2Point2Bearing{B}(x1::B) where {B <: IIF.SamplableBelief} = new{B}(x1, [P2P2BearingReuse() for i in 1:Threads.nthreads()])
end
Pose2Point2Bearing(x1::B) where {B <: IIF.SamplableBelief} = Pose2Point2Bearing{B}(x1)
function getSample(pp2br::Pose2Point2Bearing, N::Int=1)
  return (reshape(rand(pp2br.bearing, N),1,N), )
end
# define the conditional probability constraint
function (pp2br::Pose2Point2Bearing)(res::Array{Float64},
                                     userdata::FactorMetadata,
                                     idx::Int,
                                     meas::Tuple,
                                     xi::Array{Float64,2},
                                     lm::Array{Float64,2}  )
  #
  reuse = pp2br.reuse[Threads.threadid()]
  reuse.measvec[1] = cos(meas[1][idx] + xi[3,idx])
  reuse.measvec[2] = sin(meas[1][idx] + xi[3,idx])

  @simd for i in 1:2
    reuse.predvec[i] = lm[i,idx]-xi[i,idx]
  end
  reuse.resid .= reuse.measvec
  reuse.resid .-= reuse.predvec
  reuse.resid .^= 2
  res[1] = reuse.resid[1]
  res[1] += reuse.resid[2]
  return res[1]
end


# Packing and Unpacking
mutable struct PackedPose2Point2Bearing <: IncrementalInference.PackedInferenceType
    bearstr::String
    PackedPose2Point2Bearing() = new()
    PackedPose2Point2Bearing(s1::AS) where {AS <: AbstractString} = new(string(s1))
end
function convert(::Type{PackedPose2Point2Bearing}, d::Pose2Point2Bearing{B}) where {B <: IIF.SamplableBelief}
  return PackedPose2Point2Bearing(string(d.bearing))
end
function convert(::Type{Pose2Point2Bearing}, d::PackedPose2Point2Bearing)
  Pose2Point2Bearing(extractdistribution(d.bearstr))
end
