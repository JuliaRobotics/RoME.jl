##==============================================================================
## Partial Pose3 Priors
##==============================================================================
# Partial prior constraint on Z, Roll and Pitch
# ------------------------------------------------------------------------------

"""
    $(TYPEDEF)

Partial prior belief on Z, Roll, and Pitch of a `Pose3`.
"""
mutable struct PriorPose3ZRP{T1<:SamplableBelief,T2<:SamplableBelief} <: IncrementalInference.AbstractPrior
  z::T1
  rp::T2
  partial::Tuple{Int,Int,Int}
  PriorPose3ZRP{T1,T2}() where {T1, T2} = new{T1,T2}()
  PriorPose3ZRP{T1,T2}(z::T1,rp::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(z, rp, (3,4,5))
end

PriorPose3ZRP(z::T1,rp::T2) where {T1 <: SamplableBelief, T2 <: SamplableBelief} = PriorPose3ZRP{T1,T2}(z, rp)

getManifold(::PriorPose3ZRP) = ProductGroup(ProductManifold(TranslationGroup(1), SpecialOrthogonal(3)))

#FIXME update
function getSample(cfo::CalcFactor{<:PriorPose3ZRP}, N::Int=1)
  return ([[rand(cfo.factor.z); rand(cfo.factor.rp)] for _=1:N], )
end


"""
    $TYPEDEF

Serialization type of `PriorPose3ZRP`.
"""
mutable struct PackedPriorPose3ZRP <: IIF.PackedInferenceType
  zdata::String
  rpdata::String
  PackedPriorPose3ZRP() = new()
  PackedPriorPose3ZRP(x1::AS,x2::AS) where {AS <:AbstractString} = new(x1,x2)
end
function convert(::Type{PriorPose3ZRP}, d::PackedPriorPose3ZRP)
  # TODO: Change out for extractdistributionJson
  PriorPose3ZRP( convert(SamplableBelief, d.zdata), convert(SamplableBelief, d.rpdata)  )
end
function convert(::Type{PackedPriorPose3ZRP}, d::PriorPose3ZRP)
  PackedPriorPose3ZRP( convert(PackedSamplableBelief, d.z), convert(PackedSamplableBelief, d.rp) )
end


function compare(a::PriorPose3ZRP, b::PriorPose3ZRP; tol::Float64=1e-10)
  TP = true
  TP = TP && compareDensity(a.rp, b.rp)
  TP = TP && compareDensity(a.z, b.z)
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  return TP
end



##==============================================================================
## Partial Pose3 Pose3 Factors
##==============================================================================
# Partial pairwise constraint between poses X,Y,Yaw
# ------------------------------------------------------------------------------

"""
    $TYPEDEF

Partial factor between XY and Yaw of two Pose3 variables.

"""
mutable struct Pose3Pose3XYYaw{T1 <: SamplableBelief,T2 <: SamplableBelief} <: IIF.AbstractManifoldMinimize
  xy::T1
  yaw::T2
  partial::Tuple{Int,Int,Int}
  Pose3Pose3XYYaw{T1,T2}() where {T1, T2} = new()
  Pose3Pose3XYYaw{T1,T2}(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new(xy, yaw, (1,2,6))
end

Pose3Pose3XYYaw(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} =  Pose3Pose3XYYaw{T1,T2}(xy, yaw)

function getSample(cfo::CalcFactor{<:Pose3Pose3XYYaw}, N::Int=1)
  return ([[rand(cfo.factor.xy);rand(cfo.factor.yaw)] for _=1:N], )
end

getManifold(::Pose3Pose3XYYaw) = SpecialEuclidean(2)

function (cfo::CalcFactor{<:Pose3Pose3XYYaw})(meas,
                                                wXi,
                                                wXj  )
  #
  wXjhat = SE2(wXi[[1;2;6]]) * SE2(meas[1:3])
  jXjhat = SE2(wXj[[1;2;6]]) \ wXjhat
  return se2vee(jXjhat)
end


"""
    $TYPEDEF

Serialization type of Pose3Pose3XYYaw.
"""
mutable struct PackedPose3Pose3XYYaw <: IncrementalInference.PackedInferenceType
  xydata::String
  yawdata::String
  PackedPose3Pose3XYYaw() = new()
  PackedPose3Pose3XYYaw(xy::String, yaw::String) = new(xy, yaw)
end

function convert(::Type{Pose3Pose3XYYaw}, d::PackedPose3Pose3XYYaw)
  return Pose3Pose3XYYaw( convert(SamplableBelief, d.xydata), convert(SamplableBelief, d.yawdata) )
end

function convert(::Type{PackedPose3Pose3XYYaw}, d::Pose3Pose3XYYaw)
  return PackedPose3Pose3XYYaw( convert(PackedSamplableBelief, d.xy), convert(PackedSamplableBelief, d.yaw) )
end

function compare(a::Pose3Pose3XYYaw, b::Pose3Pose3XYYaw; tol::Float64=1e-10)
  TP = true
  TP = TP && compareDensity(a.xy, b.xy)
  TP = TP && compareDensity(a.yaw, b.yaw)
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  return TP
end




##==============================================================================
##
##==============================================================================
#TODO what is this, can it be removed? Moved here from above

"""
Converter: PriorPose3ZRP::Dict{String, Any} -> PriorPose3ZRP

DevNotes
- FIXME drop _evalType approach, use convert(SamplableBelief, obj) instead?
"""
function convert(::Type{RoME.PriorPose3ZRP}, fact::Dict{String, Any})
    rp = fact["measurement"][1]
    z = fact["measurement"][2]
    # FIXME drop _evalType
    rp = convert(_evalType(rp["distType"]), rp)
    z = convert(_evalType(z["distType"]), z)
    return PriorPose3ZRP(rp, z)
end

"""
Converter: PriorPose3ZRP::Dict{String, Any} -> PriorPose3ZRP
"""
function convert(::Type{Dict{String, Any}}, fact::RoME.PriorPose3ZRP)
    pf = Dict{String, Any}(
        "measurement" => [
            convert(Dict{String, Any}, fact.rp),
            convert(Dict{String, Any}, fact.z)
        ],
        "factorType" => "PriorPose3ZRP"
    )
    return pf
end

"""
    $SIGNATURES

Converter: Dict{String, Any} -> Pose3Pose3XYYaw
"""
function convert(::Type{Dict{String, Any}}, fact::RoME.Pose3Pose3XYYaw)
    pf = Dict{String, Any}(
        "measurement" => [
            convert(Dict{String, Any}, fact.xy),
            convert(Dict{String, Any}, fact.yaw)
        ],
        "factorType" => "Pose3Pose3XYYaw"
    )
    return pf
end

"""
    $SIGNATURES

Converter: Pose3Pose3XYYaw -> Dict{String, Any}

DevNotes
- FIXME stop using _evalType, see DFG #590
"""
function convert(::Type{RoME.Pose3Pose3XYYaw}, fact::Dict{String, Any})
    xy = fact["measurement"][1]
    yaw = fact["measurement"][2]
    xy = convert(_evalType(xy["distType"]), xy)
    yaw = convert(_evalType(yaw["distType"]), yaw)
    return PriorPose3ZRP(xy, yaw)
end
