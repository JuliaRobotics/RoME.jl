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

getManifold(::PriorPose3ZRP) = SpecialEuclidean(3)

#FIXME update to also only one measurement
function getSample(cf::CalcFactor{<:PriorPose3ZRP})
  
  M = getManifold(cf.factor)

  ϵ = identity_element(M)

  #Rotation part: roll and pitch
  r,p = rand(cf.factor.rp)
  R = _Rot.RotYX(p, r) #TODO confirm RotYX(p,r) or RotXY(r,p)
  
  # Translation part: Z
  T = [0; 0; rand(cf.factor.z)]

  return ProductRepr(T, R)
end


"""
    $TYPEDEF

Serialization type of `PriorPose3ZRP`.
"""
mutable struct PackedPriorPose3ZRP <: AbstractPackedFactor
  zdata::String
  rpdata::String
  # PackedPriorPose3ZRP() = new()
  # PackedPriorPose3ZRP(x1::AS,x2::AS) where {AS <:AbstractString} = new(x1,x2)
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

```
wR2 = wR1*1R2 = wR1*(1Rψ*Rθ*Rϕ)
wRz = wR1*1Rz
zRz = wRz \\ wR(Δψ)

M_R = SO(3)
δ(α,β,γ) = vee(M_R, R_0, log(M_R, R_0, zRz))

M = SE(3)
p0 = identity_element(M)
δ(x,y,z,α,β,γ) = vee(M, p0, log(M, p0, zRz))
```

"""
struct Pose3Pose3XYYaw{T <: SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T
  # partial::Tuple{Int,Int,Int,Int,Int}
  partial::Tuple{Int,Int,Int}
end
Pose3Pose3XYYaw(xy::SamplableBelief, yaw::SamplableBelief) = error("Pose3Pose3XYYaw(xy::SamplableBelief, yaw::SamplableBelief) where {T1 <: , T2 <: IIF.SamplableBelief} is deprecated, use one belief")

# Lie exponentials (pqr) are all three affected by changes in Yaw
# Pose3Pose3XYYaw(z::SamplableBelief) = Pose3Pose3XYYaw(z, (1,2,4,5,6))   # (1,2,6))
Pose3Pose3XYYaw(z::SamplableBelief) = Pose3Pose3XYYaw(z, (1,2,6))

function getSample(cf::CalcFactor{<:Pose3Pose3XYYaw})
  return sampleTangent(getManifold(cf.factor), cf.factor.Z)
end

getManifold(::Pose3Pose3XYYaw) = SpecialEuclidean(2)


## NOTE, Yaw only works if you assume a preordained global reference point, such as identity_element(Pose3)
function (cfo::CalcFactor{<:Pose3Pose3XYYaw})(X, wTp, wTq )
  #
  M = SpecialEuclidean(2)

  rx = normalize(view(wTp.parts[2],1:2, 1))
  R = SA[rx[1] -rx[2];
         rx[2]  rx[1]]
  p = ProductRepr(view(wTp.parts[1], 1:2), R)

  rx = normalize(view(wTq.parts[2],1:2, 1))
  R = SA[rx[1] -rx[2];
         rx[2]  rx[1]]
  q = ProductRepr(view(wTq.parts[1], 1:2), R)

  q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) 
  #TODO allocalte for vee! see Manifolds #412, fix for AD
  Xc = zeros(3)
  vee!(M, Xc, q, log(M, q, q̂))
  return Xc

end

## Old code and other ideas:
  # # Yaw is around the positive z axis
  # ##
  # # new Manifolds code (tests failing)
  # M3 = SpecialEuclidean(3)
  # M2 = SpecialEuclidean(2)
  
  # e3 = identity_element(M3, wTp)
  # e2 = identity_element(M2)


  # wRpψ_2  = Rotations.RotZ( TU.convert(Euler,  TU.SO3(wTp.parts[2])).Y )[1:2,1:2]
  # wTp_2  = ProductRepr(wTp.parts[1][1:2],   wRpψ_2)
  # wTqhat = Manifolds.compose(M2, wTp_2, exp(M2, e2, X))

  # wRqψ_2  = Rotations.RotZ( TU.convert(Euler,  TU.SO3(wTq.parts[2])).Y )[1:2,1:2]
  # wTq_2  = ProductRepr(wTq.parts[1][1:2],   wRqψ_2)
  # qhatTq = Manifolds.compose(M2, inv(M2, wTqhat), wTq_2)
  
  # #TODO allocate for vee! see Manifolds #412, fix for AD
  # Xc = zeros(3)
  # vee!(M2, Xc, e2, log(M2, e2, qhatTq))
  # return Xc

  # Old YPR coordinate code, < v"0.16"
  # wXjhat = SE2(wXi[[1;2;6]]) * SE2(meas[1:3])
  # jXjhat = SE2(wXj[[1;2;6]]) \ wXjhat
  # return se2vee(jXjhat)

"""
    $TYPEDEF

Serialization type of Pose3Pose3XYYaw.
"""
mutable struct PackedPose3Pose3XYYaw <: AbstractPackedFactor
  Z::String
  # PackedPose3Pose3XYYaw() = new()
  # PackedPose3Pose3XYYaw(Z::String) = new(Z)
end

function convert(::Type{<:Pose3Pose3XYYaw}, d::PackedPose3Pose3XYYaw)
  return Pose3Pose3XYYaw( convert(SamplableBelief, d.Z))
end

function convert(::Type{PackedPose3Pose3XYYaw}, d::Pose3Pose3XYYaw)
  return PackedPose3Pose3XYYaw( convert(PackedSamplableBelief, d.Z))
end

function compare(a::Pose3Pose3XYYaw, b::Pose3Pose3XYYaw; tol::Float64=1e-10)
  TP = true
  TP = TP && compareDensity(a.Z, b.Z)
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
function convert(::Type{<:PriorPose3ZRP}, fact::Dict{String, Any})
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
function convert(::Type{Dict{String, Any}}, fact::PriorPose3ZRP)
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
function convert(::Type{Dict{String, Any}}, fact::Pose3Pose3XYYaw)
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
function convert(::Type{<:Pose3Pose3XYYaw}, fact::Dict{String, Any})
    xy = fact["measurement"][1]
    yaw = fact["measurement"][2]
    xy = convert(_evalType(xy["distType"]), xy)
    yaw = convert(_evalType(yaw["distType"]), yaw)
    return PriorPose3ZRP(xy, yaw)
end
