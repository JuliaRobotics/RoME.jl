# Pose3Pose3 evaluation functions

# ------------------------------------
"""
$(TYPEDEF)

Rigid transform factor between two Pose3 compliant variables.
"""
Base.@kwdef struct Pose3Pose3{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T = MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)]))
end

getManifold(::InstanceType{Pose3Pose3}) = getManifold(Pose3) # Manifolds.SpecialEuclidean(3)

Pose3Pose3(::UniformScaling) = Pose3Pose3()

function (cf::CalcFactor{<:Pose3Pose3})(X, p::ArrayPartition{T}, q) where T
    M = getManifold(Pose3)
    q̂ = Manifolds.compose(M, p, exp(M, getPointIdentity(M), X))

    Xc::SVector{6,T} = get_coordinates(M, q, log(M, q, q̂), DefaultOrthogonalBasis())
    return Xc
end

# function (cf::CalcFactor{<:Pose3Pose3})(X, p, q)  
#   M = cf.manifold # getManifold(Pose3)
#   ϵX = exp(M, getPointIdentity(M), X)
#   q̂ = ArrayPartition(p.x[2]*ϵX.x[1] + p.x[1], p.x[2]*ϵX.x[2])
#   Xc = vee(M, q, log(M, q, q̂))
#   return Xc
# end

#TODO Serialization

"""
$(TYPEDEF)

Serialization type for `Pose3Pose3`.
"""
Base.@kwdef struct PackedPose3Pose3 <: AbstractPackedFactor
  Z::PackedSamplableBelief
end
function convert(::Type{Pose3Pose3}, packed::PackedPose3Pose3)
  return Pose3Pose3( convert(SamplableBelief, packed.Z) )
end
function convert(::Type{PackedPose3Pose3}, obj::Pose3Pose3)
  return PackedPose3Pose3( convert(PackedSamplableBelief, obj.Z) )
end




#

##
Base.@kwdef struct Pose3Pose3RotOffset{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T = MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)]))
end

getManifold(::InstanceType{Pose3Pose3RotOffset}) = getManifold(Pose3) # Manifolds.SpecialEuclidean(3)


# measurement is in frame a, for example imu frame
# p and q is in frame b, for example body frame
# bRa is the rotation to get a in the b frame 
# measurement in frame a is converted to frame b and used to calculate the error
function (cf::CalcFactor{<:Pose3Pose3RotOffset})(aX, p, q, bRa)
    M = getManifold(Pose3Pose3RotOffset)
    # measurement in frame a, input is tangent, can also use vector transport 
    a_m = exp(M, getPointIdentity(M), aX)
    b_m = ArrayPartition(a_m.x[1], bRa * a_m.x[2]) 

    q̂ = Manifolds.compose(M, p, b_m)
    return vee(M, q, log(M, q, q̂)) # coordinates
end
