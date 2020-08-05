
export VelPose2VelPose2, PackedVelPose2VelPose2


"""
$(TYPEDEF)
"""
mutable struct VelPose2VelPose2{T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} <: AbstractRelativeFactorMinimize
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



function compare(a::VelPose2VelPose2, b::VelPose2VelPose2; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && RoME.compareDensity(a.Zpose, b.Zpose)
  TP = TP && RoME.compareDensity(a.Zvel, b.Zvel)
  TP = TP && norm(a.reuseres - b.reuseres) < tol
  return TP
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
