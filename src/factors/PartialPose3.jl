# Partial Pose3 constraints

"""
    $(TYPEDEF)

Partial prior belief on Z, Roll, and Pitch of a `Pose3`.
"""
mutable struct PriorPose3ZRP{T1,T2} <: IncrementalInference.AbstractPrior where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  z::T1
  rp::T2
  partial::Tuple{Int,Int,Int}
  PriorPose3ZRP{T1,T2}() where {T1, T2} = new{T1,T2}()
  PriorPose3ZRP{T1,T2}(z::T1,rp::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(z, rp, (3,4,5))
end
PriorPose3ZRP(z::T1,rp::T2) where {T1 <: SamplableBelief, T2 <: SamplableBelief} = PriorPose3ZRP{T1,T2}(z, rp)
function getSample(pprz::PriorPose3ZRP, N::Int=1)
  return ([rand(pprz.z,N)[:]';rand(pprz.rp,N)], )
end


"""
    $SIGNATURES

Deprecated in favor of PriorPose3ZRP
"""
function convert(::Type{RoME.PriorPose3ZRP}, fact::Dict{String, Any})
    rp = fact["measurement"][1]
    z = fact["measurement"][2]
    rp = convert(_evalType(rp["distType"]), rp)
    z = convert(_evalType(z["distType"]), z)
    return PriorPose3ZRP(rp, z)
end

"""
Converter: PriorPose3ZRP::Dict{String, Any} -> PriorPose3ZRP

Deprecated in favor of PriorPose3ZRP
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
    $(TYPEDEF)

Partial prior belief on Roll Pitch and Z of a `Pose3` variable.
"""
mutable struct PartialPriorRollPitchZ{T1,T2} <: IncrementalInference.AbstractPrior where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  rp::T1
  z::T2
  partial::Tuple{Int,Int,Int}
  PartialPriorRollPitchZ{T1,T2}() where {T1, T2} = new{T1,T2}()
  PartialPriorRollPitchZ{T1,T2}(rp::T1,z::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(rp, z, (3,4,5))
end
PartialPriorRollPitchZ(rp::T1,z::T2) where {T1 <: SamplableBelief, T2 <: SamplableBelief} = PartialPriorRollPitchZ{T1,T2}(rp, z)

function getSample(cfo::CalcFactor{<:PartialPriorRollPitchZ}, N::Int=1)
  return ([rand(cfo.factor.z,N)[:]';rand(cfo.factor.rp,N)], )
end

"""
    $TYPEDEF

Serialization type of `PartialPriorRollPitchZ`.
"""
mutable struct PackedPartialPriorRollPitchZ <: IIF.PackedInferenceType
  rpdata::String
  zdata::String
  PackedPartialPriorRollPitchZ() = new()
  PackedPartialPriorRollPitchZ(x1::AS,x2::AS) where {AS <:AbstractString} = new(x1,x2)
end
function convert(::Type{PartialPriorRollPitchZ}, d::PackedPartialPriorRollPitchZ)
  # TODO: Change out for extractdistributionJson
  PartialPriorRollPitchZ( convert(SamplableBelief, d.rpdata), convert(SamplableBelief, d.zdata)  )
end
function convert(::Type{PackedPartialPriorRollPitchZ}, d::PartialPriorRollPitchZ)
  PackedPartialPriorRollPitchZ( convert(PackedSamplableBelief, d.rp), convert(PackedSamplableBelief, d.z) )
end

"""
Converter: PartialPriorRollPitchZ::Dict{String, Any} -> PartialPriorRollPitchZ

DevNotes
- FIXME drop _evalType approach, use convert(SamplableBelief, obj) instead?
"""
function convert(::Type{RoME.PartialPriorRollPitchZ}, fact::Dict{String, Any})
    rp = fact["measurement"][1]
    z = fact["measurement"][2]
    # FIXME drop _evalType
    rp = convert(_evalType(rp["distType"]), rp)
    z = convert(_evalType(z["distType"]), z)
    return PartialPriorRollPitchZ(rp, z)
end

"""
Converter: PartialPriorRollPitchZ::Dict{String, Any} -> PartialPriorRollPitchZ
"""
function convert(::Type{Dict{String, Any}}, fact::RoME.PartialPriorRollPitchZ)
    pf = Dict{String, Any}(
        "measurement" => [
            convert(Dict{String, Any}, fact.rp),
            convert(Dict{String, Any}, fact.z)
        ],
        "factorType" => "PartialPriorRollPitchZ"
    )
    return pf
end



function compare(a::PartialPriorRollPitchZ, b::PartialPriorRollPitchZ; tol::Float64=1e-10)
  TP = true
  TP = TP && compareDensity(a.rp, b.rp)
  TP = TP && compareDensity(a.z, b.z)
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  return TP
end




# Partial pairwise constraint between poses X,Y,Yaw
# ------------------------------------------------------------------------------

# TODO duplicate name Pose3Pose3XYYaw

"""
    $TYPEDEF

Partial factor between XY and Yaw of two Pose3 variables.

To be deprecated: use Pose3Pose3XYYaw instead.
"""
mutable struct PartialPose3XYYaw{T1,T2} <: AbstractRelativeMinimize where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  xy::T1
  yaw::T2
  partial::Tuple{Int,Int,Int}
  PartialPose3XYYaw{T1,T2}() where {T1, T2} = new()
  PartialPose3XYYaw{T1,T2}(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new(xy, yaw, (1,2,6))
end
PartialPose3XYYaw(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} =  PartialPose3XYYaw{T1,T2}(xy, yaw)

function getSample(cfo::CalcFactor{<:PartialPose3XYYaw}, N::Int=1)
  return ([rand(cfo.factor.xy,N);rand(cfo.factor.yaw,N)[:]'], )
end
function (cfo::CalcFactor{<:PartialPose3XYYaw})(res::AbstractVector{<:Real},
                                                meas,
                                                wXi,
                                                wXj  )
  #
  wXjhat = SE2(wXi[[1;2;6]]) * SE2(meas[1:3])
  jXjhat = SE2(wXj[[1;2;6]]) \ wXjhat
  se2vee!(res, jXjhat)
  res'*res
end

"""
    $TYPEDEF

Partial factor between XY and Yaw of two Pose3 variables.
"""
mutable struct Pose3Pose3XYYaw{T1,T2} <: AbstractRelativeMinimize where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  xy::T1
  yaw::T2
  partial::Tuple{Int,Int,Int}
  Pose3Pose3XYYaw{T1,T2}() where {T1, T2} = new()
  Pose3Pose3XYYaw{T1,T2}(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new(xy, yaw, (1,2,6))
end
Pose3Pose3XYYaw(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} =  Pose3Pose3XYYaw{T1,T2}(xy, yaw)

function getSample(pxyy::Pose3Pose3XYYaw, N::Int=1)
  return ([rand(pxyy.xy,N);rand(pxyy.yaw,N)[:]'], )
end
function (pxyy::Pose3Pose3XYYaw)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            wXi::Array{Float64,2},
            wXj::Array{Float64,2}  )
  #
  wXjhat = SE2(wXi[[1;2;6],idx]) * SE2(meas[1][1:3,idx])
  jXjhat = SE2(wXj[[1;2;6],idx]) \ wXjhat
  se2vee!(res, jXjhat)
  res'*res
end

"""
    $TYPEDEF

Serialization type of PartialPose3XYYaw.
"""
mutable struct PackedPartialPose3XYYaw <: IncrementalInference.PackedInferenceType
  xydata::String
  yawdata::String
  PackedPartialPose3XYYaw() = new()
  PackedPartialPose3XYYaw(xy::String, yaw::String) = new(xy, yaw)
end
function convert(::Type{PartialPose3XYYaw}, d::PackedPartialPose3XYYaw)
  return PartialPose3XYYaw( convert(SamplableBelief, d.xydata), convert(SamplableBelief, d.yawdata) )
end
function convert(::Type{PackedPartialPose3XYYaw}, d::PartialPose3XYYaw)
  return PackedPartialPose3XYYaw( convert(PackedSamplableBelief, d.xy), convert(PackedSamplableBelief, d.yaw) )
end

"""
    $SIGNATURES

Converter: Dict{String, Any} -> PartialPose3XYYaw
"""
function convert(::Type{Dict{String, Any}}, fact::RoME.PartialPose3XYYaw)
    pf = Dict{String, Any}(
        "measurement" => [
            convert(Dict{String, Any}, fact.xy),
            convert(Dict{String, Any}, fact.yaw)
        ],
        "factorType" => "PartialPose3XYYaw"
    )
    return pf
end

"""
    $SIGNATURES

Converter: PartialPose3XYYaw -> Dict{String, Any}

DevNotes
- FIXME stop using _evalType, see DFG #590
"""
function convert(::Type{RoME.PartialPose3XYYaw}, fact::Dict{String, Any})
    xy = fact["measurement"][1]
    yaw = fact["measurement"][2]
    xy = convert(_evalType(xy["distType"]), xy)
    yaw = convert(_evalType(yaw["distType"]), yaw)
    return PartialPriorRollPitchZ(xy, yaw)
end


function compare(a::PartialPose3XYYaw, b::PartialPose3XYYaw; tol::Float64=1e-10)
  TP = true
  TP = TP && compareDensity(a.xy, b.xy)
  TP = TP && compareDensity(a.yaw, b.yaw)
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  return TP
end














#
