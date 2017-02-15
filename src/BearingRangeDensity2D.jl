

#-------------------------------------------------------------------------------
# bearing and range available

type Pose2DPoint2DBearingRangeDensity <: IncrementalInference.FunctorPairwise
    bearing::BallTreeDensity
    range::BallTreeDensity
    Pose2DPoint2DBearingRangeDensity() = new()
    Pose2DPoint2DBearingRangeDensity(x1::BallTreeDensity,x2::BallTreeDensity) = new(x1,x2)
end
function getSample(pp2br::Pose2DPoint2DBearingRangeDensity, N::Int=1)
  b = KernelDensityEstimate.sample(pp2br.bearing, N)[1]
  r = KernelDensityEstimate.sample(pp2br.range, N)[1]
  return ([b;r], )
end
# define the conditional probability constraint
function (pp2br::Pose2DPoint2DBearingRangeDensity)(res::Array{Float64},
        idx::Int,
        meas::Tuple{Array{Float64,2}},
        xi::Array{Float64,2},
        lm::Array{Float64,2} )
  #
  # @show size(lm), size(xi), size(meas[1])
  # @show idx
  # @show meas[1][:,idx]
  # @show xi[:, idx]
  # @show lm[:,idx]
  res[1] = lm[1,idx] - (meas[1][2,idx]*cos(meas[1][1,idx]+xi[3,idx]) + xi[1,idx])
  res[2] = lm[2,idx] - (meas[1][2,idx]*sin(meas[1][1,idx]+xi[3,idx]) + xi[2,idx])
  nothing
end





# better to use bearingrange with [uniform bearing], numerical solving issue on 1D
type Pose2DPoint2DRangeDensity <: IncrementalInference.FunctorPairwise
    range::BallTreeDensity
    Pose2DPoint2DRangeDensity() = new()
    Pose2DPoint2DRangeDensity(x...) = new(x[1])
end
function (pp2r::Pose2DPoint2DRangeDensity)(res::Array{Float64},
      idx::Int,
      meas::Tuple{Array{Float64,2}, Array{Float64,1}}, # from getSample
      xi::Array{Float64,2},
      lm::Array{Float64,2}  )
  #
  # TODO -- still need to add multi-hypotheses support here

  # @show size(lm), size(xi), size(meas), size(meas[1]), size(meas[2])
  XX = lm[1,idx] - (meas[1][1,idx]*cos(meas[2][idx]+xi[3,idx]) + xi[1,idx])
  YY = lm[2,idx] - (meas[1][1,idx]*sin(meas[2][idx]+xi[3,idx]) + xi[2,idx])
  res[1] = XX^2 + YY^2
  nothing
end
function getSample(pp2r::Pose2DPoint2DRangeDensity, N::Int=1)
  return (KernelDensityEstimate.sample(pp2r.range, N)[1],  2*pi*rand(N))
end






type PackedPose2DPoint2DBearingRangeDensity <: IncrementalInference.PackedInferenceType
    bpts::Vector{Float64} # 0rotations, 1translation in each column
    bbw::Vector{Float64}
    rpts::Vector{Float64}
    rbw::Vector{Float64}
    PackedPose2DPoint2DBearingRangeDensity() = new()
    PackedPose2DPoint2DBearingRangeDensity(x...) = new(x[1], x[2], x[3], x[4])
end
function convert(::Type{Pose2DPoint2DBearingRangeDensity}, d::PackedPose2DPoint2DBearingRangeDensity)
  return Pose2DPoint2DBearingRangeDensity(
    kde!( EasyMessage(d.bpts',d.bbw)), kde!(EasyMessage(d.rpts', d.rbw) )
  )
end
function convert(::Type{PackedPose2DPoint2DBearingRangeDensity}, d::Pose2DPoint2DBearingRangeDensity)
  return PackedPose2DPoint2DBearingRangeDensity(getPoints(d.bearing)[1,:], getBW(d.bearing)[:,1],
                                                getPoints(d.range)[1,:], getBW(d.range)[:,1] )
end



type PackedPose2DPoint2DRangeDensity <: IncrementalInference.PackedInferenceType
    rpts::Vector{Float64} # 0rotations, 1translation in each column
    rbw::Vector{Float64}
    PackedPose2DPoint2DRangeDensity() = new()
    PackedPose2DPoint2DRangeDensity(x...) = new(x[1], x[2])
end
function convert(::Type{Pose2DPoint2DRangeDensity}, d::PackedPose2DPoint2DRangeDensity)
  return Pose2DPoint2DRangeDensity( kde!(EasyMessage(d.rpts', d.rbw)) )
end
function convert(::Type{PackedPose2DPoint2DRangeDensity}, d::Pose2DPoint2DRangeDensity)
  return PackedPose2DPoint2DRangeDensity( getPoints(d.range)[1,:], getBW(d.range)[:,1] )
end










#
