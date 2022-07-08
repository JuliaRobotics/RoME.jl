##==============================================================================
## Partial Pose3 Priors
##==============================================================================
# Partial prior constraint on Z, Roll and Pitch
# ------------------------------------------------------------------------------

"""
    $(TYPEDEF)

Partial prior belief on Z, Roll, and Pitch of a `Pose3`.
"""
Base.@kwdef struct PriorPose3ZRP{T1<:SamplableBelief,T2<:SamplableBelief} <: IncrementalInference.AbstractPrior
  z::T1
  rp::T2
  partial::Tuple{Int,Int,Int} = (3,4,5)
end
PriorPose3ZRP(z::SamplableBelief,rp::SamplableBelief) = PriorPose3ZRP(;z, rp)

# TODO should be dim 3 manifold
getManifold(::PriorPose3ZRP) = getManifold(Pose3) # SpecialEuclidean(3)

#FIXME update to also only one measurement
function getSample(cf::CalcFactor{<:PriorPose3ZRP})

  #Rotation part: roll and pitch
  r,p = rand(cf.factor.rp)
  R = _Rot.RotYX(p, r) #TODO confirm RotYX(p,r) or RotXY(r,p)
  
  # Translation part: Z
  T = [0; 0; rand(cf.factor.z)]

  return ArrayPartition(T, R)
end


"""
    $TYPEDEF

Serialization type of `PriorPose3ZRP`.
"""
Base.@kwdef struct PackedPriorPose3ZRP <: AbstractPackedFactor
  zdata::PackedSamplableBelief
  rpdata::PackedSamplableBelief
end
function convert(::Type{PriorPose3ZRP}, d::PackedPriorPose3ZRP)
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

getManifold(::Pose3Pose3XYYaw) = SpecialEuclidean(2)


## NOTE, Yaw only works if you assume a preordained global reference point, such as identity_element(Pose3)
function (cfo::CalcFactor{<:Pose3Pose3XYYaw})(X, wTp, wTq )
  #
  M = SpecialEuclidean(2)

  rx = normalize(view(wTp.x[2],1:2, 1))
  R = SA[rx[1] -rx[2];
         rx[2]  rx[1]]
  p = ArrayPartition(view(wTp.x[1], 1:2), R)

  rx = normalize(view(wTq.x[2],1:2, 1))
  R = SA[rx[1] -rx[2];
         rx[2]  rx[1]]
  q = ArrayPartition(view(wTq.x[1], 1:2), R)

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
Base.@kwdef struct PackedPose3Pose3XYYaw <: AbstractPackedFactor
  Z::PackedSamplableBelief
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




#
