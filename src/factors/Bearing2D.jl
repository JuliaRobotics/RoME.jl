
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
struct Pose2Point2Bearing{B <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
    bearing::B
    reuse::Vector{P2P2BearingReuse}
    Pose2Point2Bearing{B}() where B = new{B}()
    Pose2Point2Bearing{B}(x1::B) where {B <: IIF.SamplableBelief} = new{B}(x1, [P2P2BearingReuse() for i in 1:Threads.nthreads()])
end
Pose2Point2Bearing(x1::B) where {B <: IIF.SamplableBelief} = Pose2Point2Bearing{B}(x1)

function getSample(cfo::CalcFactor{<:Pose2Point2Bearing}, N::Int=1)
  return (reshape(rand(cfo.factor.bearing, N),1,N), )
end
# define the conditional probability constraint
function (cfo::CalcFactor{<:Pose2Point2Bearing})( res::AbstractVector{<:Real},
                                                  meas,
                                                  xi,
                                                  lm  )
  #
  reuse = cfo.factor.reuse[Threads.threadid()]
  reuse.measvec[1] = cos(meas[1] + xi[3])
  reuse.measvec[2] = sin(meas[1] + xi[3])

  # @simd for i in 1:2
    reuse.predvec[1:2] .= lm[1:2].-xi[1:2]
  # end
  normalize!(reuse.predvec)
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
  return PackedPose2Point2Bearing(convert(PackedSamplableBelief, d.bearing))
end
function convert(::Type{Pose2Point2Bearing}, d::PackedPose2Point2Bearing)
  Pose2Point2Bearing(convert(SamplableBelief, d.bearstr))
end
