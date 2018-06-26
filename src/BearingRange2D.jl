# Bearing and Range constraints for 2D


# better to use bearingrange with [uniform bearing], numerical solving issue on 1D
mutable struct Pose2DPoint2DRange <: IncrementalInference.FunctorPairwise
    Zij::Vector{Float64} # bearing and range hypotheses as columns
    Cov::Float64
    W::Vector{Float64}
    Pose2DPoint2DRange() = new()
    Pose2DPoint2DRange(x...) = new(x[1],x[2],x[3])
end
function getSample(pp2::Pose2DPoint2DRange, N::Int=1)
  return (pp2.Cov*randn(1,N),  2*pi*rand(N))
end
function (pp2r::Pose2DPoint2DRange)(res::Array{Float64},
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

mutable struct Pose2DPoint2DBearingRange{B <: Distributions.Distribution, R <: Distributions.Distribution} <: IncrementalInference.FunctorPairwise
    # Zij::Array{Float64,2} # bearing and range hypotheses as columns
    # Cov::Array{Float64,2}
    # W::Array{Float64,1}
    bearing::B
    range::R
    Pose2DPoint2DBearingRange{B,R}() where {B,R} = new{B,R}()
    Pose2DPoint2DBearingRange(x1::B,x2::R) where {B,R} = new{B,R}(x1,x2)
    # Pose2DPoint2DBearingRange{B,R}(x1::B,x2::R) where {B,R} = new{B,R}(x1,x2)
end
function getSample(pp2br::Pose2DPoint2DBearingRange, N::Int=1)
  b = rand(pp2br.bearing, N)
  r = rand(pp2br.range, N)
  return ([b';r'], )
end
# define the conditional probability constraint
function (pp2br::Pose2DPoint2DBearingRange)(res::Array{Float64},
        userdata,
        idx::Int,
        meas::Tuple{Array{Float64,2}},
        xi::Array{Float64,2},
        lm::Array{Float64,2} )
  #
  res[1] = lm[1,idx] - (meas[1][2,idx]*cos(meas[1][1,idx]+xi[3,idx]) + xi[1,idx])
  res[2] = lm[2,idx] - (meas[1][2,idx]*sin(meas[1][1,idx]+xi[3,idx]) + xi[2,idx])
  nothing
end



# Support for database based solving

passTypeThrough(d::FunctionNodeData{Pose2DPoint2DRange}) = d

mutable struct PackedPose2DPoint2DBearingRange <: IncrementalInference.PackedInferenceType
    # bmu::Float64 # 0rotations, 1translation in each column
    # bsig::Float64
    # rmu::Float64
    # rsig::Float64
    bearstr::String
    rangstr::String
    PackedPose2DPoint2DBearingRange() = new()
    # PackedPose2DPoint2DBearingRange(x1, x2, x3, x4) = new(x1, x2, x3, x4)
    PackedPose2DPoint2DBearingRange(s1::AS, s2::AS) where {AS <: AbstractString} = new(string(s1),string(s2))
end


function convert(::Type{PackedPose2DPoint2DBearingRange}, d::Pose2DPoint2DBearingRange{Normal{T}, Normal{T}}) where T
  return PackedPose2DPoint2DBearingRange(string(d.bearing), string(d.range))
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types
function convert(::Type{Pose2DPoint2DBearingRange}, d::PackedPose2DPoint2DBearingRange)
  Pose2DPoint2DBearingRange(extractdistribution(d.bearstr), extractdistribution(d.rangstr))
end


mutable struct Pose2DPoint2DBearingRangeMH{B <: Distributions.Distribution, R <: Distributions.Distribution} <: IncrementalInference.FunctorPairwise
    bearing::B
    range::R
    hypothesis::Distributions.Categorical
    Pose2DPoint2DBearingRangeMH{B,R}() where {B,R} = new{B,R}()
    Pose2DPoint2DBearingRangeMH(x1::B,x2::R, w::Distributions.Categorical) where {B,R} = new{B,R}(x1,x2,w)
    Pose2DPoint2DBearingRangeMH(x1::B,x2::R, w::Vector{Float64}=Float64[1.0;]) where {B,R} = new{B,R}(x1,x2,Categorical(w))
end
function getSample(pp2br::Pose2DPoint2DBearingRangeMH, N::Int=1)::Tuple{Array{Float64,2}, Vector{Int}}
  b = rand(pp2br.bearing, N)
  r = rand(pp2br.range, N)
  s = rand(pp2br.hypothesis, N)
  return ([b';r'], s)
end
# define the conditional probability constraint
function (pp2br::Pose2DPoint2DBearingRangeMH)(res::Array{Float64},
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

mutable struct PackedPose2DPoint2DBearingRangeMH <: IncrementalInference.PackedInferenceType
    bearstr::String
    rangstr::String
    hypostr::String
    PackedPose2DPoint2DBearingRangeMH() = new()
    PackedPose2DPoint2DBearingRangeMH(s1::AS, s2::AS, s3::AS) where {AS <: AbstractString} = new(string(s1),string(s2),string(s3))
end
function convert(::Type{PackedPose2DPoint2DBearingRangeMH}, d::Pose2DPoint2DBearingRangeMH{Normal{T}, Normal{T}}) where T
  return PackedPose2DPoint2DBearingRangeMH(string(d.bearing), string(d.range), string(d.hypothesis))
end
# TODO -- should not be resorting to string, consider specialized code for parametric distribution types
function convert(::Type{Pose2DPoint2DBearingRangeMH}, d::PackedPose2DPoint2DBearingRangeMH)
  Pose2DPoint2DBearingRangeMH(extractdistribution(d.bearstr), extractdistribution(d.rangstr), extractdistribution(d.hypostr))
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















# ------------------------------------------------------
