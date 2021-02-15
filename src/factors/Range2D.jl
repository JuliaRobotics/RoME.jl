


"""
$(TYPEDEF)
"""
mutable struct Point2Point2Range{D <: IIF.SamplableBelief} <: IncrementalInference.AbstractRelativeMinimize
  Z::D
  Point2Point2Range{D}() where {D} = new{D}()
  Point2Point2Range{D}(d::D) where {D <: IIF.SamplableBelief} = new{D}(d)
end
Point2Point2Range(d::D) where {D <: IIF.SamplableBelief} = Point2Point2Range{D}(d)
function getSample(cfo::CalcFactor{<:Point2Point2Range}, N::Int=1)
  return (reshape(rand(cfo.factor.Z,N),1,N),  2*pi*rand(N))
end
function (cfo::CalcFactor{<:Point2Point2Range})(rho, theta, xi, lm)
  #
  XX = lm[1] - (rho[1]*cos(theta[1]) + xi[1])
  YY = lm[2] - (rho[1]*sin(theta[1]) + xi[2])
  return XX^2 + YY^2
end

#Parametric only function fro Point2Point2Range 
function (cfo::CalcFactor{<:Point2Point2Range})(rho, xi, lm)
  # Basically `EuclidDistance`
  return [rho[1] - norm(lm .- xi)]
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
  return PackedPoint2Point2Range(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Point2Point2Range}, d::PackedPoint2Point2Range)
  return Point2Point2Range(convert(SamplableBelief, d.str))
end



"""
    $TYPEDEF

Range only measurement from Pose2 to Point2 variable.
"""
mutable struct Pose2Point2Range{T} <: IIF.AbstractRelativeMinimize
  Z::T
  partial::Tuple{Int,Int}
  Pose2Point2Range{T}() where T = new()
  Pose2Point2Range{T}(Z::T) where {T <: IIF.SamplableBelief} = new{T}(Z, (1,2))
end
Pose2Point2Range(Z::T) where {T <: IIF.SamplableBelief} = Pose2Point2Range{T}(Z)

function getSample(cfo::CalcFactor{<:Pose2Point2Range}, N::Int=1)
  return (reshape(rand(cfo.factor.Z,N),1,N) ,  2*pi*rand(N))
end
function (pp2r::CalcFactor{<:Pose2Point2Range})(rho,
                                                theta,
                                                xi,
                                                lm  )
  #
  # DONE in IIF -- still need to add multi-hypotheses support here
  # this is the noisy range
  XX = lm[1] - (rho[1]*cos(theta[1]) + xi[1])
  YY = lm[2] - (rho[1]*sin(theta[1]) + xi[2])
  return sqrt(XX^2 + YY^2)  
end



mutable struct PackedPose2Point2Range  <: IncrementalInference.PackedInferenceType
  str::String
  PackedPose2Point2Range() = new()
  PackedPose2Point2Range(s::AS) where {AS <: AbstractString} = new(s)
end
function convert(::Type{PackedPose2Point2Range}, d::Pose2Point2Range)
  return PackedPose2Point2Range(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{Pose2Point2Range}, d::PackedPose2Point2Range)
  return Pose2Point2Range(convert(SamplableBelief, d.str))
end
