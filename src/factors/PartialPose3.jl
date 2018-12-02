# Partial Pose3 constraints

"""
$(TYPEDEF)
"""
mutable struct PriorPose3ZRP{T1,T2} <: IncrementalInference.FunctorSingleton where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  z::T1
  rp::T2
  partial::Tuple
  PriorPose3ZRP{T1,T2}() where {T1, T2} = new{T1,T2}()
  PriorPose3ZRP{T1,T2}(z::T1,rp::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(z, rp, (3,4,5))
end
PriorPose3ZRP(z::T1,rp::T2) where {T1 <: SamplableBelief, T2 <: SamplableBelief} = PriorPose3ZRP{T1,T2}(z, rp)
function getSample(pprz::PriorPose3ZRP, N::Int=1)
  return ([rand(pprz.z,N)[:]';rand(pprz.rp,N)], )
end


"""
Converter: PartialPriorRollPitchZ::Dict{String, Any} -> PartialPriorRollPitchZ

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
"""
mutable struct PartialPriorRollPitchZ{T1,T2} <: IncrementalInference.FunctorSingleton where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  rp::T1
  z::T2
  partial::Tuple
  PartialPriorRollPitchZ{T1,T2}() where {T1, T2} = new{T1,T2}()
  PartialPriorRollPitchZ{T1,T2}(rp::T1,z::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(rp, z, (3,4,5))
end
PartialPriorRollPitchZ(rp::T1,z::T2) where {T1 <: SamplableBelief, T2 <: SamplableBelief} = PartialPriorRollPitchZ{T1,T2}(rp, z)
function getSample(pprz::PartialPriorRollPitchZ, N::Int=1)
  return ([rand(pprz.z,N)[:]';rand(pprz.rp,N)], )
end

mutable struct PackedPartialPriorRollPitchZ <: IncrementalInference.PackedInferenceType
  rpdata::String
  zdata::String
  PackedPartialPriorRollPitchZ() = new()
  PackedPartialPriorRollPitchZ(x1::AS,x2::AS) where {AS <:AbstractString} = new(x1,x2)
end
function convert(::Type{PartialPriorRollPitchZ}, d::PackedPartialPriorRollPitchZ)
  # TODO: Change out for extractdistributionJson
  PartialPriorRollPitchZ( extractdistribution(d.rpdata), extractdistribution(d.zdata)  )
end
function convert(::Type{PackedPartialPriorRollPitchZ}, d::PartialPriorRollPitchZ)
  PackedPartialPriorRollPitchZ( string(d.rp), string(d.z) )
end

"""
Converter: PartialPriorRollPitchZ::Dict{String, Any} -> PartialPriorRollPitchZ
"""
function convert(::Type{RoME.PartialPriorRollPitchZ}, fact::Dict{String, Any})
    rp = fact["measurement"][1]
    z = fact["measurement"][2]
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


function compare(a::Normal, b::Normal; tol::Float64=1e-10)
  TP = true
  TP = TP && abs(a.μ - b.μ) < tol
  TP = TP && abs(a.σ - b.σ) < tol
  TP
end

function compare(a::PartialPriorRollPitchZ, b::PartialPriorRollPitchZ; tol::Float64=1e-10)
  TP = true
  TP = TP && compare(a.rp, b.rp)
  TP = TP && compare(a.z, b.z)
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  return TP
end




# Partial pairwise constraint between poses X,Y,Yaw
# ------------------------------------------------------------------------------

# TODO duplicate name Pose3Pose3XYYaw

mutable struct PartialPose3XYYaw{T1,T2} <: FunctorPairwise where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  xy::T1
  yaw::T2
  partial::Tuple
  PartialPose3XYYaw{T1,T2}() where {T1, T2} = new()
  PartialPose3XYYaw{T1,T2}(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new(xy, yaw, (1,2,6))
end
PartialPose3XYYaw(xy::T1, yaw::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} =  PartialPose3XYYaw{T1,T2}(xy, yaw)

function getSample(pxyy::PartialPose3XYYaw, N::Int=1)
  return ([rand(pxyy.xy,N);rand(pxyy.yaw,N)[:]'], )
end
function (pxyy::PartialPose3XYYaw)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            wXi::Array{Float64,2},
            wXj::Array{Float64,2}  )
  #
  wXjhat = SE2(wXi[[1;2;6],idx]) * SE2(meas[1][1:3,idx])
  jXjhat = SE2(wXj[[1;2;6],idx]) \ wXjhat
  se2vee!(res, jXjhat)
  nothing
end

mutable struct Pose3Pose3XYYaw{T1,T2} <: FunctorPairwise where {T1 <: SamplableBelief, T2 <: SamplableBelief}
  xy::T1
  yaw::T2
  partial::Tuple
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
  nothing
end


mutable struct PackedPartialPose3XYYaw <: IncrementalInference.PackedInferenceType
  xydata::String
  yawdata::String
  PackedPartialPose3XYYaw() = new()
  PackedPartialPose3XYYaw(xy::String, yaw::String) = new(xy, yaw)
end
function convert(::Type{PartialPose3XYYaw}, d::PackedPartialPose3XYYaw)
  return PartialPose3XYYaw( extractdistribution(d.xydata), extractdistribution(d.yawdata) )
end
function convert(::Type{PackedPartialPose3XYYaw}, d::PartialPose3XYYaw)
  return PackedPartialPose3XYYaw( string(d.xy), string(d.yaw) )
end

"""
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
Converter: PartialPose3XYYaw -> Dict{String, Any}
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
  TP = TP && compare(a.xy, b.xy)
  TP = TP && compare(a.yaw, b.yaw)
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  return TP
end










mutable struct PartialPose3XYYawNH <: IncrementalInference.FunctorPairwiseNH
  xyy::Distributions.MvNormal
  partial::Tuple
  nullhypothesis::Distributions.Categorical
  PartialPose3XYYawNH() = new()
  PartialPose3XYYawNH(xyy::MvNormal, vh::Vector{Float64}) = new(xyy, (1,2,6),  Distributions.Categorical(vh))
end
function getSample(pxyy::PartialPose3XYYawNH, N::Int=1)
  return (rand(pxyy.xyy,N), )
end
function (pxyy::PartialPose3XYYawNH)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            wXi::Array{Float64,2},
            wXj::Array{Float64,2}  )
  #
  wXjhat = SE2(wXi[[1;2;6],idx])*SE2(meas[1][:,idx]) #*SE2(pp2.Zij[:,1])*SE2(meas[1][:,idx])
  jXjhat = SE2(wXj[[1;2;6],idx]) \ wXjhat
  se2vee!(res, jXjhat)
  nothing
end



mutable struct PackedPartialPose3XYYawNH <: IncrementalInference.PackedInferenceType
  vecZij::Array{Float64,1} # 3translations, 3rotation
  vecCov::Array{Float64,1}
  nullhypothesis::Vector{Float64}
  PackedPartialPose3XYYawNH() = new()
  PackedPartialPose3XYYawNH(x1::Vector{Float64}, x2::Array{Float64}, x3::Vector{Float64}) = new(x1, x2[:], x3)
end
function convert(::Type{PartialPose3XYYawNH}, d::PackedPartialPose3XYYawNH)
  return PartialPose3XYYawNH( Distributions.MvNormal(d.vecZij,
               reshapeVec2Mat(d.vecCov, 3)), d.nullhypothesis  )
end
function convert(::Type{PackedPartialPose3XYYawNH}, d::PartialPose3XYYawNH)
  return PackedPartialPose3XYYawNH(d.xyy.μ, d.xyy.Σ.mat, d.nullhypothesis.p )
end


function compare(a::PartialPose3XYYawNH, b::PartialPose3XYYawNH; tol::Float64=1e-10)
  TP = true
  TP = TP && norm(a.xyy.μ-b.xyy.μ) < tol
  TP = TP && norm(a.xyy.Σ.mat[:]-b.xyy.Σ.mat[:]) < tol
  TP = TP && norm(collect(a.partial)-collect(b.partial)) < tol
  TP = TP && norm(a.nullhypothesis.p-b.nullhypothesis.p) < tol
  return TP
end






#
