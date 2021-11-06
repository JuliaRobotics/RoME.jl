
export SpecialEuclidean


"""
$(TYPEDEF)

Rigid transform between two Pose2's, assuming (x,y,theta).

DevNotes
- Maybe with Manifolds.jl, `{T <: IIF.SamplableBelief, S, R, P}`

Related

[`Pose3Pose3`](@ref), [`Point2Point2`](@ref), [`MutablePose2Pose2Gaussian`](@ref), [`DynPose2`](@ref), [`InertialPose3`](ref)
"""
struct Pose2Pose2{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Z::T
end
# convenience and default constructor
Pose2Pose2() = Pose2Pose2(MvNormal(Diagonal([1.0; 1.0; 1.0])))
 
DFG.getManifold(::Pose2Pose2) = Manifolds.SpecialEuclidean(2)

Pose2Pose2(::UniformScaling) = Pose2Pose2() # MvNormal(zeros(3),LinearAlgebra.diagm([1.0;1.0;1.0])) )

#TODO remove to use default
function getSample(cf::CalcFactor{<:Pose2Pose2}) 
  
  M = getManifold(cf.factor)
  ϵ = getPointIdentity(Pose2)

  X = sampleTangent(M, cf.factor.Z, ϵ)
  return X
end


# Assumes X is a tangent vector
function (cf::CalcFactor{<:Pose2Pose2})(X, p, q)
    @assert X isa ProductRepr "Pose2Pose2 expects measurement sample X to be a Manifolds tangent vector, not coordinate or point representation.  Got X=$X"
    M = getManifold(Pose2)
    q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
    #TODO allocalte for vee! see Manifolds #412, fix for AD
    Xc = zeros(3)
    vee!(M, Xc, q, log(M, q, q̂))
    return Xc
end
  


# NOTE, serialization support -- will be reduced to macro in future
# ------------------------------------

"""
$(TYPEDEF)
"""
mutable struct PackedPose2Pose2  <: IIF.PackedInferenceType
  datastr::String
  # PackedPose2Pose2() = new()
  # PackedPose2Pose2(x::String) = new(x)
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
  return Pose2Pose2(convert(SamplableBelief, d.datastr))
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
  return PackedPose2Pose2(convert(PackedSamplableBelief, d.Z))
end



# FIXME, rather have separate compareDensity functions
function compare(a::Pose2Pose2,b::Pose2Pose2; tol::Float64=1e-10)
  return compareDensity(a.Z, b.Z)
  # TP = true
  # TP = TP && norm(a.z.μ-b.z.μ) < (tol + 1e-5)
  # TP = TP && norm(a.z.Σ.mat-b.z.Σ.mat) < tol
  # return TP
end

## Deprecated
# function Pose2Pose2(mean::Array{Float64,1}, cov::Array{Float64,2}, w::Vector{Float64})
#   @warn "Pose2Pose2(mu,cov,w) is deprecated in favor of Pose2Pose2(T(...)) -- use for example Pose2Pose2(MvNormal(mu, cov))"
#   Pose2Pose2(MvNormal(mean, cov))
# end
