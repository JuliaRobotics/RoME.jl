# Pose2 functions for Robot Motion Estimate

struct Pose2 <: IncrementalInference.InferenceVariable
  dims::Int
  labels::Vector{String}
  Pose2() = new(3, String["POSE";])
end

# TODO - move to IncrementalInference
struct Prior{T} <: IncrementalInference.FunctorSingleton where {T <: Distribution}
  z::T
end
getSample(s::Prior, N::Int=1) = (rand(s.z,N), )


mutable struct Pose2Pose2{T} <: IncrementalInference.FunctorPairwise where {T <: Distributions.Distribution}
  z::T
  Pose2Pose2{T}() where {T <: Distribution} = new{T}()
  Pose2Pose2(z1::T) where {T <: Distribution} = new{T}(z1)
end
Pose2Pose2(mean::Array{Float64,1}, cov::Array{Float64,2}) = Pose2Pose2(MvNormal(mean, cov))
Pose2Pose2(mean::Array{Float64,1}, cov::Array{Float64,2}, w::Vector{Float64}) = Pose2Pose2(MvNormal(mean, cov))
getSample(s::Pose2Pose2, N::Int=1) = (rand(s.z,N), )
function (s::Pose2Pose2)(res::Array{Float64},
            userdata,
            idx::Int,
            meas::Tuple,
            wxi::Array{Float64,2},
            wxj::Array{Float64,2}  )
  #
  wXjhat = SE2(wxi[:,idx])*SE2(meas[1][:,idx])
  jXjhat = SE2(wxj[:,idx]) \ wXjhat
  se2vee!(res, jXjhat)
  nothing
end


mutable struct PriorPose2 <: IncrementalInference.FunctorSingleton
    Zi::Array{Float64,2}
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    PriorPose2() = new()
    PriorPose2(x...) = new(x[1], x[2], x[3])
end
function getSample(p2::PriorPose2, N::Int=1)
  return (rand(MvNormal(p2.Zi[:,1],p2.Cov),N), )
end



function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  TP = true
  TP = TP && norm(a.Zi-b.Zi) < tol
  TP = TP && norm(a.Cov-b.Cov) < tol
  TP = TP && norm(a.W-b.W) < tol
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






# NOTE, for database support -- will be reduced to macro in future
# ------------------------------------


mutable struct PackedPriorPose2  <: IncrementalInference.PackedInferenceType
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int
    vecCov::Array{Float64,1}
    dimc::Int
    W::Array{Float64,1}
    PackedPriorPose2() = new()
    PackedPriorPose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  Zi = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPose2(Zi, Cov, d.W)
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  v1 = d.Zi[:];
  v2 = d.Cov[:];
  return PackedPriorPose2(v1,size(d.Zi,1),
                          v2,size(d.Cov,1),
                          d.W)
end



# --------------------------------------------





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
