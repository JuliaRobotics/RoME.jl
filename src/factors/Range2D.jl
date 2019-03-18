


"""
$(TYPEDEF)
"""
mutable struct Point2Point2Range{D <: IIF.SamplableBelief} <: IncrementalInference.FunctorPairwiseMinimize
  Z::D
  Point2Point2Range{D}() where {D} = new{D}()
  Point2Point2Range{D}(d::D) where {D <: IIF.SamplableBelief} = new{D}(d)
end
Point2Point2Range(d::D) where {D <: IIF.SamplableBelief} = Point2Point2Range{D}(d)
function getSample(pp2::Point2Point2Range{T}, N::Int=1) where {T <: IIF.SamplableBelief}
  return (reshape(rand(pp2.Z,N),1,N),  2*pi*rand(N))
end
function (pp2r::Point2Point2Range{T})(
            res::Array{Float64},
            userdata::FactorMetadata,
            idx::Int,
            meas::Tuple,
            xi::Array{Float64,2},
            lm::Array{Float64,2} ) where {T <: IIF.SamplableBelief}
  #
  z = meas[1][1,idx]
  XX = lm[1,idx] - (z*cos(meas[2][idx]) + xi[1,idx])
  YY = lm[2,idx] - (z*sin(meas[2][idx]) + xi[2,idx])
  res[1] = XX^2 + YY^2
  # @show "fnc", pointer(res), res
  res[1]
end
# import RoME: Point2Point2Range



"""
    $TYPEDEF

Range only measurement from Pose2 to Point2 variable.
"""
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





passTypeThrough(d::FunctionNodeData{Point2Point2Range}) = d

"""
$(TYPEDEF)
"""
mutable struct PackedPoint2Point2Range  <: IncrementalInference.PackedInferenceType
  str::String
  PackedPoint2Point2Range() = new()
  PackedPoint2Point2Range(s::AS) where {AS <: AbstractString} = new(s)
end
function convert(::Type{PackedPoint2Point2Range}, d::Point2Point2Range)
  return PackedPoint2Point2Range(string(d.Z))
end
function convert(::Type{Point2Point2Range}, d::PackedPoint2Point2Range)
  return Point2Point2Range(extractdistribution(d.str))
end
