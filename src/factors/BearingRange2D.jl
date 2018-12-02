# Bearing and Range constraints for 2D

# better to use bearingrange with [uniform bearing], numerical solving issue on 1D
mutable struct Pose2Point2Range{T} <: IncrementalInference.FunctorPairwise
  Z::T
  Pose2Point2Range{T}() where T = new()
  Pose2Point2Range{T}(Z::T) where {T <: IIF.SamplableBelief} = new{T}(Z)
end
Pose2Point2Range(Z::T) where {T <: IIF.SamplableBelief} = Pose2Point2Range{T}(Z)

function getSample(pp2::Pose2Point2Range, N::Int=1)
  return (reshape(rand(pp2.Z,N),1,N) ,  2*pi*rand(N))
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
  z = meas[1][1,idx]
  XX = lm[1,idx] - (z*cos(meas[2][idx]) + xi[1,idx])
  YY = lm[2,idx] - (z*sin(meas[2][idx]) + xi[2,idx])
  res[1] = XX^2 + YY^2
  nothing
end

#-------------------------------------------------------------------------------
# bearing and range available


mutable struct Pose2Point2BearingRange{B <: IIF.SamplableBelief, R <: IIF.SamplableBelief} <: IncrementalInference.FunctorPairwise
    bearing::B
    range::R
    Pose2Point2BearingRange{B,R}() where {B,R} = new{B,R}()
    Pose2Point2BearingRange{B,R}(x1::B,x2::R) where {B <: IIF.SamplableBelief,R <: IIF.SamplableBelief} = new{B,R}(x1,x2)
end
Pose2Point2BearingRange(x1::B,x2::R) where {B <: IIF.SamplableBelief,R <: IIF.SamplableBelief} = Pose2Point2BearingRange{B,R}(x1,x2)
function getSample(pp2br::Pose2Point2BearingRange, N::Int=1)
  smpls = zeros(2, N)
  smpls[1,:] = rand(pp2br.bearing, N)[:]
  smpls[2,:] = rand(pp2br.range, N)[:]
  return (smpls,)
end
# define the conditional probability constraint
function (pp2br::Pose2Point2BearingRange)(res::Array{Float64},
        userdata::FactorMetadata,
        idx::Int,
        meas::Tuple{Array{Float64,2}},
        xi::Array{Float64,2},
        lm::Array{Float64,2} )
  #
  res[1] = lm[1,idx] - (meas[1][2,idx]*cos(meas[1][1,idx]+xi[3,idx]) + xi[1,idx])
  res[2] = lm[2,idx] - (meas[1][2,idx]*sin(meas[1][1,idx]+xi[3,idx]) + xi[2,idx])
  nothing
end

# import RoME: Pose2Point2BearingRange


# Support for database based solving

passTypeThrough(d::FunctionNodeData{Pose2Point2Range}) = d

mutable struct PackedPose2Point2BearingRange <: IncrementalInference.PackedInferenceType
    bearstr::String
    rangstr::String
    PackedPose2Point2BearingRange() = new()
    PackedPose2Point2BearingRange(s1::AS, s2::AS) where {AS <: AbstractString} = new(string(s1),string(s2))
end

function convert(::Type{PackedPose2Point2BearingRange}, d::Pose2Point2BearingRange{B, R}) where {B <: IIF.SamplableBelief, R <: IIF.SamplableBelief}
  return PackedPose2Point2BearingRange(string(d.bearing), string(d.range))
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert(::Type{Pose2Point2BearingRange}, d::PackedPose2Point2BearingRange)
 # where {B <: IIF.SamplableBelief, R <: IIF.SamplableBelief}
  Pose2Point2BearingRange(extractdistribution(d.bearstr), extractdistribution(d.rangstr))
end







#-------------------------------------------------------------------------------
# bearing only available

# this factor type is still a work in progress
mutable struct Pose2Point2Bearing{B <: Distributions.Distribution} <: IncrementalInference.FunctorPairwiseMinimize
    bearing::B
    Pose2Point2Bearing{B}() where B = new{B}()
    Pose2Point2Bearing{B}(x1::B) where {B <: Distributions.Distribution} = new{B}(x1)
end
Pose2Point2Bearing(x1::B) where {B <: Distributions.Distribution} = Pose2Point2Bearing{B}(x1)
function getSample(pp2br::Pose2Point2Bearing, N::Int=1)
  return (reshape(rand(pp2br.bearing, N),1,N), )
end
# define the conditional probability constraint
function (pp2br::Pose2Point2Bearing)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple,
            xi::Array{Float64,2},
            lm::Array{Float64,2}  )
  #
  res[1] = ( meas[1][idx] - atan(lm[2,idx]-xi[2,idx], lm[1,idx]-xi[1,idx]) )^2
  return res[1]
end










# ------------------------------------------------------
