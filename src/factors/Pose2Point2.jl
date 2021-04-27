
export Pose2Point2, PackedPose2Point2

#-------------------------------------------------------------------------------
#

"""
    $TYPEDEF

Bearing and Range constraint from a Pose2 to Point2 variable.
"""
struct Pose2Point2{T <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
    Zij::T
    partial::Tuple{Int,Int}
    # empty constructor
    # Pose2Point2{T}() where {T} = new{T}()
    # # regular constructor
    # Pose2Point2{T}(x1::T) where {T <: IIF.SamplableBelief} = new{T}(x1)
end
# convenience and default constructor
Pose2Point2(x1::T=MvNormal(zeros(2),LinearAlgebra.diagm([0.01;0.01]))) where {T <: IIF.SamplableBelief} = Pose2Point2{T}(x1, (1,2))

# prescribed sampling function
function getSample(cfo::CalcFactor{<:Pose2Point2}, N::Int=1)
  smpls = zeros(2, N)
  smpls[1:2,:] .= rand(cfo.factor.Zij, N)

  return (smpls,)
end

function IIF.getParametricMeasurement(s::Pose2Point2{<:MvNormal})

  return IIF.getParametricMeasurement(s.Zij)

end

# define the conditional probability constraint
function (cfo::CalcFactor{<:Pose2Point2})(meas,
                                          wXi,
                                          wLj )
  #

  wLj_pred = SE2(wXi)*SE2([meas;0.0])
  return  wLj .- se2vee(wLj_pred)[1:2]
end


## Serialization support

mutable struct PackedPose2Point2 <: IncrementalInference.PackedInferenceType
    Zij::String
    PackedPose2Point2() = new()
    PackedPose2Point2(s1::AS) where {AS <: AbstractString} = new(string(s1))
end

function convert(::Type{PackedPose2Point2}, obj::Pose2Point2{T}) where {T <: IIF.SamplableBelief}
  return PackedPose2Point2(convert(PackedSamplableBelief, obj.Zij))
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert(::Type{Pose2Point2}, packed::PackedPose2Point2)
  Pose2Point2(convert(SamplableBelief, packed.Zij))
end
