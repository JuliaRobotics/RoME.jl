# Bearing and Range constraints for 2D

SamplableBelief = Union{Distributions.Distribution, KernelDensityEstimate.BallTreeDensity}

# function samplebelief(d::Distribution, N::Int=1)
#   rand(d, N)
# end
# function samplebelief(d::BallTreeDensity, N::Int=1)
#   KernelDensityEstimate.sample(d, N)[1]
# end

# better to use bearingrange with [uniform bearing], numerical solving issue on 1D
mutable struct Pose2Point2Range{T} <: IncrementalInference.FunctorPairwise
  Z::T
  Pose2Point2Range{T}() = new()
  Pose2Point2Range{T}(Z::T) where {T <: SamplableBelief} = new{T}(Z)
end
Pose2Point2Range(Z::T) where {T <: SamplableBelief} = Pose2Point2Range{T}(Z)

function Pose2DPoint2DRange(x1::Vector{T},x2::Array{T,2},x3) where {T <: Real}
  warn("Pose2Point2Range(mu,cov,w) is being deprecated in favor of Pose2Point2Range(T(...)), such as Pose2Point2Range(MvNormal(mu, cov))")
  Pose2Point2Range(MvNormal(x1, x2))
end

function getSample(pp2::Pose2Point2Range, N::Int=1)
  return (pp2.Cov*randn(1,N),  2*pi*rand(N))
end
function (pp2r::Pose2Point2Range)(res::Array{Float64},
                                    userdata,
                                    idx::Int,
                                    meas::Tuple{Array{Float64,2}, Array{Float64,1}}, # from getSample
                                    xi::Array{Float64,2},
                                    lm::Array{Float64,2}  )
  #
  # DONE in IIF -- still need to add multi-hypotheses support here
  # this is the noisy range
  z = pp2r.Zij[1]+meas[1][1,idx]
  XX = lm[1,idx] - (z*cos(meas[2][idx]) + xi[1,idx])
  YY = lm[2,idx] - (z*sin(meas[2][idx]) + xi[2,idx])
  res[1] = XX^2 + YY^2
  nothing
end

#-------------------------------------------------------------------------------
# bearing and range available



mutable struct Pose2Point2BearingRange{B <: SamplableBelief, R <: SamplableBelief} <: IncrementalInference.FunctorPairwise
    bearing::B
    range::R
    Pose2Point2BearingRange{B,R}() where {B,R} = new{B,R}()
    Pose2Point2BearingRange{B,R}(x1::B,x2::R) where {B <: SamplableBelief,R <: SamplableBelief} = new{B,R}(x1,x2)
end
Pose2Point2BearingRange(x1::B,x2::R) where {B <: SamplableBelief,R <: SamplableBelief} = Pose2Point2BearingRange{B,R}(x1,x2)
function getSample(pp2br::Pose2Point2BearingRange, N::Int=1)
  b = rand(pp2br.bearing, N)
  r = rand(pp2br.range, N)
  return (b, r)
end
# define the conditional probability constraint
function (pp2br::Pose2Point2BearingRange)(res::Array{Float64},
        userdata::FactorMetadata,
        idx::Int,
        meas::Tuple{Vector{Float64}, Vector{Float64}},
        xi::Array{Float64,2},
        lm::Array{Float64,2} )
  #
  res[1] = lm[1,idx] - (meas[2][idx]*cos(meas[1][idx]+xi[3,idx]) + xi[1,idx])
  res[2] = lm[2,idx] - (meas[2][idx]*sin(meas[1][idx]+xi[3,idx]) + xi[2,idx])
  nothing
end



# Support for database based solving

passTypeThrough(d::FunctionNodeData{Pose2Point2Range}) = d

mutable struct PackedPose2Point2BearingRange <: IncrementalInference.PackedInferenceType
    bearstr::String
    rangstr::String
    PackedPose2Point2BearingRange() = new()
    PackedPose2Point2BearingRange(s1::AS, s2::AS) where {AS <: AbstractString} = new(string(s1),string(s2))
end

function convert(::Type{PackedPose2Point2BearingRange}, d::Pose2Point2BearingRange{B, R}) where {B <: SamplableBelief, R <: SamplableBelief}
  return PackedPose2Point2BearingRange(string(d.bearing), string(d.range))
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert(::Type{Pose2Point2BearingRange}, d::PackedPose2Point2BearingRange)
 # where {B <: SamplableBelief, R <: SamplableBelief}
  Pose2Point2BearingRange(extractdistribution(d.bearstr), extractdistribution(d.rangstr))
end


mutable struct Pose2Point2BearingRangeMH{B <: Distributions.Distribution, R <: Distributions.Distribution} <: IncrementalInference.FunctorPairwise
    bearing::B
    range::R
    hypothesis::Distributions.Categorical
    Pose2Point2BearingRangeMH{B,R}() where {B,R} = new{B,R}()
    Pose2Point2BearingRangeMH(x1::B,x2::R, w::Distributions.Categorical) where {B,R} = new{B,R}(x1,x2,w)
    Pose2Point2BearingRangeMH(x1::B,x2::R, w::Vector{Float64}=Float64[1.0;]) where {B,R} = new{B,R}(x1,x2,Categorical(w))
end
function getSample(pp2br::Pose2Point2BearingRangeMH, N::Int=1)::Tuple{Array{Float64,2}, Vector{Int}}
  b = rand(pp2br.bearing, N)
  r = rand(pp2br.range, N)
  s = rand(pp2br.hypothesis, N)
  return ([b';r'], s)
end
# define the conditional probability constraint
function (pp2br::Pose2Point2BearingRangeMH)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple{Array{Float64,2}, Vector{Int}},
            xi::Array{Float64,2},
            lms... )::Void  # ::Array{Float64,2}
  #
  warn("Older interface, not analytically correct.")
  res[1] = lms[meas[2][idx]][1,idx] - (meas[1][2,idx]*cos(meas[1][1,idx]+xi[3,idx]) + xi[1,idx])
  res[2] = lms[meas[2][idx]][2,idx] - (meas[1][2,idx]*sin(meas[1][1,idx]+xi[3,idx]) + xi[2,idx])
  nothing
end

mutable struct PackedPose2Point2BearingRangeMH <: IncrementalInference.PackedInferenceType
    bearstr::String
    rangstr::String
    hypostr::String
    PackedPose2Point2BearingRangeMH() = new()
    PackedPose2Point2BearingRangeMH(s1::AS, s2::AS, s3::AS) where {AS <: AbstractString} = new(string(s1),string(s2),string(s3))
end
function convert(::Type{PackedPose2Point2BearingRangeMH}, d::Pose2Point2BearingRangeMH{Normal{T}, Normal{T}}) where T
  return PackedPose2Point2BearingRangeMH(string(d.bearing), string(d.range), string(d.hypothesis))
end
# TODO -- should not be resorting to string, consider specialized code for parametric distribution types
function convert(::Type{Pose2Point2BearingRangeMH}, d::PackedPose2Point2BearingRangeMH)
  Pose2Point2BearingRangeMH(extractdistribution(d.bearstr), extractdistribution(d.rangstr), extractdistribution(d.hypostr))
end





#-------------------------------------------------------------------------------
# bearing only available

# this factor type is still a work in progress
mutable struct Pose2DPoint2DBearing{B <: Distributions.Distribution} <: IncrementalInference.FunctorPairwise
    bearing::B
    Pose2DPoint2DBearing{B}() where {B} = new{B}()
    Pose2DPoint2DBearing(x1::B) where {B} = new{B}(x1)
    Pose2DPoint2DBearing{B}(x1::B) where {B} = new{B}(x1)
end
function getSample(pp2br::Pose2DPoint2DBearing, N::Int=1)
  return (rand(pp2br.bearing, N), )
end
# define the conditional probability constraint
function (pp2br::Pose2DPoint2DBearing)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple,
            xi::Array{Float64,2},
            lm::Array{Float64,2}  )
  #
  res[1] = meas[1][idx] - atan2(lm[2,idx]-xi[2,idx], lm[1,idx]-xi[1,idx])
  nothing
end




# mutable struct PackedPose2Point2BearingRangeDensity <: IncrementalInference.PackedInferenceType
#     bpts::Vector{Float64} # 0rotations, 1translation in each column
#     bbw::Vector{Float64}
#     rpts::Vector{Float64}
#     rbw::Vector{Float64}
#     PackedPose2Point2BearingRangeDensity() = new()
#     PackedPose2Point2BearingRangeDensity(x...) = new(x[1], x[2], x[3], x[4])
# end
# function convert(::Type{Pose2Point2BearingRangeDensity}, d::PackedPose2Point2BearingRangeDensity)
#   return Pose2Point2BearingRangeDensity(
#     kde!( EasyMessage(d.bpts',d.bbw)), kde!(EasyMessage(d.rpts', d.rbw) )
#   )
# end
# function convert(::Type{PackedPose2Point2BearingRangeDensity}, d::Pose2Point2BearingRangeDensity)
#   return PackedPose2Point2BearingRangeDensity(getPoints(d.bearing)[1,:], getBW(d.bearing)[:,1],
#                                                 getPoints(d.range)[1,:], getBW(d.range)[:,1] )
# end










# ------------------------------------------------------
