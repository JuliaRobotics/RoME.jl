
"""
$(TYPEDEF)
"""
struct Point2 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  Point2() = new(2, String[])
end

"""
$(TYPEDEF)
"""
mutable struct PriorPoint2{T} <: IncrementalInference.FunctorSingleton where {T <: Distributions.Distribution}
  Z::T
  # W::Array{Float64,1} # TODO, deprecate the weight parameter
  PriorPoint2{T}() where T = new()
  PriorPoint2{T}(dist::T) where {T <: Distributions.Distribution} = new{T}(dist)
end
PriorPoint2(z::T) where {T <: Distributions.Distribution} = PriorPoint2{T}(z)
function PriorPoint2D(mu, cov, W)
  @warn "PriorPoint2D(mu, cov, W) is deprecated, use PriorPoint{T}(T(...)) instead -- e.g. PriorPoint2{MvNormal}(MvNormal(...) or any other Distributions.Distribution type instead."
  PriorPoint2{MvNormal{Float64}}(MvNormal(mu, cov))
end
function getSample(p2::PriorPoint2, N::Int=1)
  return (rand(p2.Z, N),)
end

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
$(TYPEDEF)
"""
mutable struct Point2Point2{D <: IIF.SamplableBelief} <: FunctorPairwise #BetweenPoses
    Zij::D
    Point2Point2{T}() where T = new{T}()
    Point2Point2{T}(x::T) where {T <: Sampleable} = new{T}(x)
end
Point2Point2(x::T) where {T <: Sampleable} = Point2Point2{T}(x)
function getSample(pp2::Point2Point2, N::Int=1)
  return (rand(pp2.Zij,N),  )
end
function (pp2r::Point2Point2{T})(
            res::Array{Float64},
            userdata::FactorMetadata,
            idx::Int,
            meas::Tuple,
            xi::Array{Float64,2},
            xj::Array{Float64,2} ) where T
  #
  res[1]  = meas[1][1,idx] - (xj[1,idx] - xi[1,idx])
  res[2]  = meas[1][2,idx] - (xj[2,idx] - xi[2,idx])
  nothing
end



"""
$(TYPEDEF)
"""
mutable struct Point2Point2WorldBearing{T} <: IncrementalInference.FunctorPairwise where {T <: IIF.SamplableBelief}
    Z::T
    rangemodel::Rayleigh
    # zDim::Tuple{Int, Int}
    Point2Point2WorldBearing{T}() where T = new{T}()
    Point2Point2WorldBearing{T}(x::T) where {T <: IIF.SamplableBelief} = new{T}(x, Rayleigh(100))
end
Point2Point2WorldBearing(x::T) where {T <: IIF.SamplableBelief} = Point2Point2WorldBearing{T}(x)
function getSample(pp2::Point2Point2WorldBearing, N::Int=1)
  sp = Array{Float64,2}(undef, 2,N)
  sp[1,:] = rand(pp2.Z,N)
  sp[2,:] = rand(pp2.rangemodel,N)
  return (sp, )
end
function (pp2r::Point2Point2WorldBearing)(
          res::Array{Float64},
          userdata::FactorMetadata,
          idx::Int,
          meas::Tuple,
          pi::Array{Float64,2},
          pj::Array{Float64,2} )
  #
  # noisy bearing measurement
  z, r = meas[1][1,idx], meas[1][2,idx]
  dx, dy = pj[1,idx]-pi[1,idx], pj[2,idx]-pi[2,idx]
  res[1] = z - atan(dy,dx)
  res[2] = r - norm([dx; dy])
  nothing
end



"""
$(TYPEDEF)
"""
mutable struct PriorPoint2DensityNH <: IncrementalInference.FunctorSingletonNH
  belief::BallTreeDensity
  nullhypothesis::Distributions.Categorical
  PriorPoint2DensityNH() = new()
  PriorPoint2DensityNH(belief, p::Distributions.Categorical) = new(belief, p)
  PriorPoint2DensityNH(belief, p::Vector{Float64}) = new(belief, Distributions.Categorical(p))
end
function getSample(p2::PriorPoint2DensityNH, N::Int=1)
  return (rand(p2.belief, N), )
end
"""
$(TYPEDEF)
"""
mutable struct PackedPriorPoint2DensityNH <: IncrementalInference.PackedInferenceType
    rpts::Vector{Float64} # 0rotations, 1translation in each column
    rbw::Vector{Float64}
    dims::Int
    nh::Vector{Float64}
    PackedPriorPoint2DensityNH() = new()
    PackedPriorPoint2DensityNH(x1,x2,x3, x4) = new(x1, x2, x3, x4)
end
function convert(::Type{PriorPoint2DensityNH}, d::PackedPriorPoint2DensityNH)
  return PriorPoint2DensityNH(
            kde!(EasyMessage( reshapeVec2Mat(d.rpts, d.dims), d.rbw)),
            Distributions.Categorical(d.nh)  )
end
function convert(::Type{PackedPriorPoint2DensityNH}, d::PriorPoint2DensityNH)
  return PackedPriorPoint2DensityNH( getPoints(d.belief)[:], getBW(d.belief)[:,1], Ndim(d.belief), d.nullhypothesis.p )
end










# ---------------------------------------------------------



"""
$(TYPEDEF)
"""
mutable struct PackedPriorPoint2  <: IncrementalInference.PackedInferenceType
    str::String
    PackedPriorPoint2() = new()
    PackedPriorPoint2(x::String) = new(x)
end

passTypeThrough(d::FunctionNodeData{Point2Point2Range}) = d

function convert(::Type{PriorPoint2}, d::PackedPriorPoint2)
  # Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  distr = extractdistribution(d.str)
  return PriorPoint2{typeof(distr)}(distr)
end
function convert(::Type{PackedPriorPoint2}, d::PriorPoint2)
  # v2 = d.mv.Σ.mat[:];
  return PackedPriorPoint2(string(d.Z))
end




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





"""
$(TYPEDEF)
"""
mutable struct PackedPoint2Point2WorldBearing  <: IncrementalInference.PackedInferenceType
    str::String
    # NOTE Not storing rangemodel which may cause inconsistencies if the implementation parameters change
    PackedPoint2Point2WorldBearing() = new()
    PackedPoint2Point2WorldBearing(x::String) = new(x)
end
function convert(::Type{PackedPoint2Point2WorldBearing}, d::Point2Point2WorldBearing)
  return PackedPoint2Point2WorldBearing( string(d.Z) )
end
function convert(::Type{Point2Point2WorldBearing}, d::PackedPoint2Point2WorldBearing)
  return Point2Point2WorldBearing( extractdistribution(d.str) )
end





"""
$(TYPEDEF)
"""
mutable struct PackedPoint2Point2 <: IncrementalInference.PackedInferenceType
    # mu::Vector{Float64}
    # sigma::Vector{Float64}
    # sdim::Int
    str::String
    PackedPoint2Point2() = new()
    PackedPoint2Point2(s::AS) where {AS <: AbstractString} = new(s)
end
function convert(::Type{Point2Point2}, d::PackedPoint2Point2)
  return Point2Point2( extractdistribution(d.str) )
end
function convert(::Type{PackedPoint2Point2}, d::Point2Point2)
  return PackedPoint2Point2( string(d.Zij) )
end
