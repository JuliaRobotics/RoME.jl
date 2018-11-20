# 2D SLAM with velocity states



"""
$(TYPEDEF)
"""
mutable struct DynPose2VelocityPrior{T1,T2} <: IncrementalInference.FunctorSingleton where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief}
  Zpose::T1
  Zvel::T2
  DynPose2VelocityPrior{T1,T2}() where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = new{T1,T2}()
  DynPose2VelocityPrior{T1,T2}(z1::T1,z2::T2) where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = new{T1,T2}(z1,z2)
end
DynPose2VelocityPrior(z1::T1,z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = DynPose2VelocityPrior{T1,T2}(z1,z2)
getSample(dp2v::DynPose2VelocityPrior{T1,T2}, N::Int=1) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = ([rand(dp2v.Zpose,N);rand(dp2v.Zvel,N)], )




# x, y, dx/dt, dy/dt
"""
$(TYPEDEF)
"""
mutable struct DynPose2 <: IncrementalInference.InferenceVariable
  ut::Int64 # microsecond time
  dims::Int
  labels::Vector{String}
  DynPose2(;ut::Int64=-9999999999, labels::Vector{<:AbstractString}=String["POSE";]) = new(ut, 5, labels)
end




"""
$(TYPEDEF)
"""
mutable struct VelPose2VelPose2{T1,T2} <: IncrementalInference.FunctorPairwiseMinimize where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief}
  Zpose::Pose2Pose2{T1} #Zpose::T1
  Zvel::T2
  reuseres::Vector{Vector{Float64}}
  VelPose2VelPose2{T1,T2}() where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}()
  VelPose2VelPose2{T1,T2}(z1::T1, z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(Pose2Pose2(z1),z2,[zeros(3) for i in 1:Threads.nthreads()])
end
VelPose2VelPose2(z1::T1, z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = VelPose2VelPose2{T1,T2}(z1, z2)
getSample(vp2vp2::VelPose2VelPose2, N::Int=1) = ([rand(vp2vp2.Zpose.z,N);rand(vp2vp2.Zvel,N)], )
function (vp2vp2::VelPose2VelPose2{T1,T2})(
                res::Array{Float64},
                userdata::FactorMetadata,
                idx::Int,
                meas::Tuple,
                Xi::Array{Float64,2},
                Xj::Array{Float64,2}  ) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief}
  #
  z = meas[1][:,idx]
  wxi, wxj = Xi[:,idx], Xj[:,idx]
  # @show z, wxi, wxj
  dt = (userdata.variableuserdata[2].ut - userdata.variableuserdata[1].ut)*1e-6   # roughly the intended use of userdata
  # fill!(vp2vp2.reuseres[Threads.threadid()], 0.0)
  vp2vp2.Zpose(vp2vp2.reuseres[Threads.threadid()], userdata, idx, meas, Xi, Xj)
  wDXij = (wxj[4:5]-wxi[4:5])
  bDXij = TransformUtils.R(-wxi[3])*wDXij
  # calculate the residual
  res[1] = sum((vp2vp2.reuseres[Threads.threadid()]).^2)
  res[1] += sum((z[4:5] - bDXij).^2)
  dx = se2vee(SE2(wxi[1:3]) \ SE2(wxj[1:3]))
  res[1] += sum((dx[1:2]/dt - 0.5*(wxj[4:5]+wxi[4:5])).^2)  # first order integration
  res[1]
end

# import RoME: VelPose2VelPose2, IIF.SamplableBelief
# import IncrementalInference: FactorMetadata





"""
$(TYPEDEF)
"""
mutable struct DynPose2Pose2{T} <: IncrementalInference.FunctorPairwise where {T <: IIF.SamplableBelief}
  Zpose::Pose2Pose2{T} #Zpose::T1
  # reuseres::Vector{Float64}
  partial::Tuple
  DynPose2Pose2{T}() where {T <: IIF.SamplableBelief} = new{T}()
  DynPose2Pose2{T}(z1::T) where {T <: IIF.SamplableBelief} = new{T}(Pose2Pose2(z1), (1,2,3))
end
DynPose2Pose2(z1::T) where {T <: IIF.SamplableBelief} = DynPose2Pose2{T}(z1)
getSample(vp2vp2::DynPose2Pose2, N::Int=1) = (rand(vp2vp2.Zpose.z,N), )
function (vp2vp2::DynPose2Pose2{T})(
                res::Array{Float64},
                userdata::FactorMetadata,
                idx::Int,
                meas::Tuple,
                Xi::Array{Float64,2},
                Xj::Array{Float64,2}  ) where {T <: IIF.SamplableBelief}
  #
  vp2vp2.Zpose(res, userdata, idx, meas, Xi, Xj)
  nothing
end




## Compare functions

function compare(a::MvNormal, b::MvNormal; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && norm(a.μ - b.μ)<tol
  TP = TP && sum(norm.(a.Σ.mat - b.Σ.mat))<tol
  return TP
end


function compare(a::DynPose2VelocityPrior, b::DynPose2VelocityPrior)::Bool
  RoME.compare(a.Zpose, b.Zpose) && RoME.compare(a.Zvel, b.Zvel)
end


function compare(a::VelPose2VelPose2, b::VelPose2VelPose2; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && RoME.compare(a.Zpose, b.Zpose)
  TP = TP && RoME.compare(a.Zvel, b.Zvel)
  TP = TP && norm(a.reuseres - b.reuseres) < tol
  return TP
end


## Packing types


"""
$(TYPEDEF)
"""
mutable struct PackedDynPose2VelocityPrior <: IncrementalInference.PackedInferenceType
  strpose::AbstractString
  strvel::AbstractString
  PackedDynPose2VelocityPrior() = new()
  PackedDynPose2VelocityPrior(z1::AS, z2::AS) where {AS <: AbstractString} = new(z1, z2)
end

function convert(::Type{PackedDynPose2VelocityPrior}, d::DynPose2VelocityPrior)
  return PackedDynPose2VelocityPrior(string(d.Zpose),string(d.Zvel))
end
function convert(::Type{DynPose2VelocityPrior}, d::PackedDynPose2VelocityPrior)
  posedistr = extractdistribution(d.strpose)
  veldistr = extractdistribution(d.strvel)
  return DynPose2VelocityPrior(posedistr, veldistr)
end



"""
$(TYPEDEF)
"""
mutable struct PackedVelPose2VelPose2 <: IncrementalInference.PackedInferenceType
  strpose::AbstractString
  strvel::AbstractString
  PackedVelPose2VelPose2() = new()
  PackedVelPose2VelPose2(z1::AS, z2::AS) where {AS <: AbstractString} = new(z1, z2)
end

function convert(::Type{PackedVelPose2VelPose2}, d::VelPose2VelPose2)
  return PackedVelPose2VelPose2(string(d.Zpose.z),string(d.Zvel))
end
function convert(::Type{VelPose2VelPose2}, d::PackedVelPose2VelPose2)
  posedistr = extractdistribution(d.strpose)
  veldistr = extractdistribution(d.strvel)
  return VelPose2VelPose2(posedistr, veldistr)
end




"""
$(TYPEDEF)
"""
mutable struct PackedDynPose2Pose2 <: IncrementalInference.PackedInferenceType
  strpose::AbstractString
  PackedDynPose2Pose2() = new()
  PackedDynPose2Pose2(z1::AS) where {AS <: AbstractString} = new(z1)
end

function convert(::Type{PackedDynPose2Pose2}, d::DynPose2Pose2)
  return PackedDynPose2Pose2(string(d.Zpose.z))
end
function convert(::Type{DynPose2Pose2}, d::PackedDynPose2Pose2)
  posedistr = extractdistribution(d.strpose)
  return DynPose2Pose2(posedistr)
end






#
