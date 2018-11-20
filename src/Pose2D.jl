
# TODO -- temporary overwriting of this function, and should be removed from here when IIF 4.0.0 is tagged.
import IncrementalInference: compare

function compare(a::IncrementalInference.GenericFunctionNodeData{T1,S},b::IncrementalInference.GenericFunctionNodeData{T2,S}) where {T1, T2, S}
  # TODO -- beef up this comparison to include the gwp
  TP = true
  TP = TP && a.fncargvID == b.fncargvID
  TP = TP && a.eliminated == b.eliminated
  TP = TP && a.potentialused == b.potentialused
  TP = TP && a.edgeIDs == b.edgeIDs
  TP = TP && a.frommodule == b.frommodule
  # TP = TP && typeof(a.fnc) == typeof(b.fnc)
  return TP
end



# Pose2 functions for Robot Motion Estimate

"""
$(TYPEDEF)
"""
struct Pose2 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  Pose2() = new(3, String["POSE";])
end

# # Done - move to IncrementalInference
# struct Prior{T} <: IncrementalInference.FunctorSingleton where {T <: Distribution}
#   z::T
# end
# getSample(s::Prior, N::Int=1) = (rand(s.z,N), )


"""
$(TYPEDEF)
"""
mutable struct Pose2Pose2{T} <: IncrementalInference.FunctorPairwise where {T <: IIF.SamplableBelief}
  z::T
  Pose2Pose2{T}() where {T <: IIF.SamplableBelief} = new{T}()
  Pose2Pose2{T}(z1::T) where {T <: IIF.SamplableBelief} = new{T}(z1)
end
Pose2Pose2(z::T) where {T <: IIF.SamplableBelief} = Pose2Pose2{T}(z)
function Pose2Pose2(mean::Array{Float64,1}, cov::Array{Float64,2})
  @warn "Pose2Pose2(mu,cov) is deprecated in favor of Pose2Pose2(T(...)) -- use for example Pose2Pose2(MvNormal(mu, cov))"
  Pose2Pose2(MvNormal(mean, cov))
end
function Pose2Pose2(mean::Array{Float64,1}, cov::Array{Float64,2}, w::Vector{Float64})
  @warn "Pose2Pose2(mu,cov,w) is deprecated in favor of Pose2Pose2(T(...)) -- use for example Pose2Pose2(MvNormal(mu, cov))"
  Pose2Pose2(MvNormal(mean, cov))
end
getSample(s::Pose2Pose2{<:IIF.SamplableBelief}, N::Int=1) = (rand(s.z,N), )
function (s::Pose2Pose2{<:IIF.SamplableBelief})(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple,
            wxi::Array{Float64,2},
            wxj::Array{Float64,2}  )
  #
  wXjhat = SE2(wxi[1:3,idx])*SE2(meas[1][1:3,idx])
  jXjhat = SE2(wxj[1:3,idx]) \ wXjhat
  se2vee!(res, jXjhat)
  nothing
end


"""
$(TYPEDEF)
"""
mutable struct PriorPose2{T} <: IncrementalInference.FunctorSingleton  where {T <: Distributions.Distribution}
    Z::T
    PriorPose2{T}() where T = new{T}()
    PriorPose2{T}(x::T) where {T <: Distributions.Distribution}  = new{T}(x)
end
PriorPose2(x::T) where {T <: Distributions.Distribution} = PriorPose2{T}(x)
function PriorPose2(mu::Array{Float64}, cov::Array{Float64,2}, W::Vector{Float64})
  @warn "PriorPose2(mu,cov,W) is deprecated in favor of PriorPose2(T(...)) -- use for example PriorPose2(MvNormal(mu, cov))"
  PriorPose2(MvNormal(mu[:], cov))
end
function getSample(p2::PriorPose2, N::Int=1)
  return (rand(p2.Z,N), )
end



function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  TP = true
  TP = TP && norm(a.Z.μ-b.Z.μ) < tol
  TP = TP && norm(a.Z.Σ.mat-b.Z.Σ.mat) < tol
  return TP
end
function compare(a::Pose2Pose2,b::Pose2Pose2; tol::Float64=1e-10)
  TP = true
  TP = TP && norm(a.z.μ-b.z.μ) < (tol + 1e-5)
  TP = TP && norm(a.z.Σ.mat-b.z.Σ.mat) < tol
  return TP
end




# Project all particles (columns) Xval with Z, that is for all  SE3(Xval[:,i])*Z
function projectParticles(Xval::Array{Float64,2}, Z::Array{Float64,2}, Cov::Array{Float64,2})
  # TODO optimize convert SE2 to a type

  r,c = size(Xval)
  RES = zeros(r,c) #*cz

  # ent, x = SE3(0), SE3(0)
  j=1
  # for j in 1:cz
  ENT = rand( MvNormal(Z[:,1], Cov), c )
    for i in 1:c
      x = SE2(Xval[1:3,i])
      dx = SE2(ENT[1:3,i])
      RES[1:r,i*j] = se2vee(x*dx)
    end
  # end
  #
  return RES
end

⊕(Xpts::Array{Float64,2}, z::Pose2Pose2) = projectParticles(Xpts, z.Zij, z.Cov)
⊕(Xvert::Graphs.ExVertex, z::Pose2Pose2) = ⊕(getVal(Xvert), z)




"""
$(TYPEDEF)
"""
mutable struct PartialPriorYawPose2{T} <: IncrementalInference.FunctorSingleton  where {T <: IIF.SamplableBelief}
    Z::T
    partial::Tuple
    PartialPriorYawPose2{T}() where T = new{T}()
    PartialPriorYawPose2{T}(x::T) where {T <: IIF.SamplableBelief}  = new{T}(x, (3,))
end
PartialPriorYawPose2(x::T) where {T <: IIF.SamplableBelief} = PartialPriorYawPose2{T}(x)

function getSample(p2::PartialPriorYawPose2, N::Int=1)
  return (reshape(rand(p2.Z,N),1,N), )
end



"""
$(TYPEDEF)
"""
mutable struct PackedPartialPriorYawPose2 <: IncrementalInference.PackedInferenceType
    Z::String
    PackedPartialPriorYawPose2() = new()
    PackedPartialPriorYawPose2(x::T) where {T <: AbstractString}  = new(x)
end

function convert(::Type{PackedPartialPriorYawPose2}, d::PartialPriorYawPose2)
  PackedPartialPriorYawPose2(string(d.Z))
end
function convert(::Type{PartialPriorYawPose2}, d::PackedPartialPriorYawPose2)
  PartialPriorYawPose2(extractdistribution(d.Z))
end




# NOTE, for database support -- will be reduced to macro in future
# ------------------------------------


"""
$(TYPEDEF)
"""
mutable struct PackedPriorPose2  <: IncrementalInference.PackedInferenceType
    # vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    # dimz::Int
    # vecCov::Array{Float64,1}
    # dimc::Int
    # W::Array{Float64,1}
    str::String
    PackedPriorPose2() = new()
    PackedPriorPose2(x::String) = new(x)
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  # v1 = d.Zi[:];
  # v2 = d.Cov[:];
  return PackedPriorPose2(string(d.Z))
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  # Zi = reshapeVec2Mat(d.vecZij,d.dimz)
  # Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  distr = extractdistribution(d.str)
  return PriorPose2{typeof(distr)}(distr)
end



# --------------------------------------------





"""
$(TYPEDEF)
"""
mutable struct PackedPose2Pose2  <: IncrementalInference.PackedInferenceType
  datastr::String
  PackedPose2Pose2() = new()
  PackedPose2Pose2(x::String) = new(x)
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
  return Pose2Pose2(extractdistribution(d.datastr))
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
  return PackedPose2Pose2(string(d.z))
end





# --------------------------------------------
