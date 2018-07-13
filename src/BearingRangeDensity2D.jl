

#-------------------------------------------------------------------------------
# bearing and range available

mutable struct Pose2Point2BearingRangeDensity <: IncrementalInference.FunctorPairwise
    bearing::BallTreeDensity
    range::BallTreeDensity
    Pose2Point2BearingRangeDensity() = new()
    Pose2Point2BearingRangeDensity(x1::BallTreeDensity,x2::BallTreeDensity) = new(x1,x2)
end
function getSample(pp2br::Pose2Point2BearingRangeDensity, N::Int=1)
  b = KernelDensityEstimate.sample(pp2br.bearing, N)[1]
  r = KernelDensityEstimate.sample(pp2br.range, N)[1]
  return ([b;r], )
end
# define the conditional probability constraint
function (pp2br::Pose2Point2BearingRangeDensity)(res::Array{Float64},
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




# better to use bearingrange with [uniform bearing], numerical solving issue on 1D
mutable struct Pose2Point2RangeDensity <: IncrementalInference.FunctorPairwiseMinimize
    range::BallTreeDensity
    Pose2Point2RangeDensity() = new()
    Pose2Point2RangeDensity(x...) = new(x[1])
end
function (pp2r::Pose2Point2RangeDensity)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple{Array{Float64,2}, Array{Float64,1}}, # from getSample
            xi::Array{Float64,2},
            lm::Array{Float64,2}  )
  # DONE in IIF -- still need to add multi-hypotheses support here
  XX = lm[1,idx] - (meas[1][1,idx]*cos(meas[2][idx]+xi[3,idx]) + xi[1,idx])
  YY = lm[2,idx] - (meas[1][1,idx]*sin(meas[2][idx]+xi[3,idx]) + xi[2,idx])
  res[1] = XX^2 + YY^2
  res[1]
end
function getSample(pp2r::Pose2Point2RangeDensity, N::Int=1)
  return (KernelDensityEstimate.sample(pp2r.range, N)[1],  2*pi*rand(N))
end






mutable struct PackedPose2Point2BearingRangeDensity <: IncrementalInference.PackedInferenceType
    bpts::Vector{Float64} # 0rotations, 1translation in each column
    bbw::Vector{Float64}
    rpts::Vector{Float64}
    rbw::Vector{Float64}
    PackedPose2Point2BearingRangeDensity() = new()
    PackedPose2Point2BearingRangeDensity(x...) = new(x[1], x[2], x[3], x[4])
end
function convert(::Type{Pose2Point2BearingRangeDensity}, d::PackedPose2Point2BearingRangeDensity)
  return Pose2Point2BearingRangeDensity(
    kde!( EasyMessage(d.bpts',d.bbw)), kde!(EasyMessage(d.rpts', d.rbw) )
  )
end
function convert(::Type{PackedPose2Point2BearingRangeDensity}, d::Pose2Point2BearingRangeDensity)
  return PackedPose2Point2BearingRangeDensity(getPoints(d.bearing)[1,:], getBW(d.bearing)[:,1],
                                                getPoints(d.range)[1,:], getBW(d.range)[:,1] )
end



mutable struct PackedPose2Point2RangeDensity <: IncrementalInference.PackedInferenceType
    rpts::Vector{Float64} # 0rotations, 1translation in each column
    rbw::Vector{Float64}
    PackedPose2Point2RangeDensity() = new()
    PackedPose2Point2RangeDensity(x...) = new(x[1], x[2])
end
function convert(::Type{Pose2Point2RangeDensity}, d::PackedPose2Point2RangeDensity)
  return Pose2Point2RangeDensity( kde!(EasyMessage(d.rpts', d.rbw)) )
end
function convert(::Type{PackedPose2Point2RangeDensity}, d::Pose2Point2RangeDensity)
  return PackedPose2Point2RangeDensity( getPoints(d.range)[1,:], getBW(d.range)[:,1] )
end










#
