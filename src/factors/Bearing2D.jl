
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
struct Pose2Point2Bearing{B <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
    bearing::B
    reuse::Vector{P2P2BearingReuse}
    Pose2Point2Bearing{B}() where B = new{B}()
    Pose2Point2Bearing{B}(x1::B) where {B <: IIF.SamplableBelief} = new{B}(x1, [P2P2BearingReuse() for i in 1:Threads.nthreads()])
end
Pose2Point2Bearing(x1::B) where {B <: IIF.SamplableBelief} = Pose2Point2Bearing{B}(x1)

getManifold(::Pose2Point2Bearing) = TranslationGroup(1)

function getSample(cfo::CalcFactor{<:Pose2Point2Bearing}, N::Int=1)
  return (rand(cfo.factor.bearing, N), )
end

function (cfo::CalcFactor{<:Pose2Point2Bearing})(Xc, p, l)
  #embed l in SE2
  M = SpecialEuclidean(2)
  q = ProductRepr(l, identity_element(SpecialOrthogonal(2), p.parts[2]))
  x,y = Manifolds.compose(M, inv(M, p), q).parts[1]
  # TODO Xc is a coordinate (ie angle), maybe change to X ϵ so2 
  # m̂ = exp(so{N}(hat(SpecialOrthogonal(N), SO{N}()[], atan(y, x))))
  # distance(m, m̂)/sqrt(2)
  return Xc[1] .- atan(y, x)
end
# define the conditional probability constraint
# function (cfo::CalcFactor{<:Pose2Point2Bearing})(meas, _xi, lm)
#   #
#   M = SpecialEuclidean(2)
#   ϵ = identity_element(M)
#   xi = vee(M, ϵ, log(M, ϵ, _xi))

#   reuse = cfo.factor.reuse[Threads.threadid()]
#   reuse.measvec[1] = cos(meas[1] + xi[3])
#   reuse.measvec[2] = sin(meas[1] + xi[3])

#   # @simd for i in 1:2
#     reuse.predvec[1:2] .= lm[1:2].-xi[1:2]
#   # end
#   normalize!(reuse.predvec)
#   reuse.resid .= reuse.measvec
#   reuse.resid .-= reuse.predvec

#   # must return result of length zDim==1 in this case
#   return sum(abs.(reuse.resid))
#   # reuse.resid .^= 2
#   # res[1] = reuse.resid[1]
#   # res[1] += reuse.resid[2]
  
# end


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
