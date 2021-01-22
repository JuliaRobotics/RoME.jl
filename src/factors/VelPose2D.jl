
export VelPose2VelPose2, PackedVelPose2VelPose2


"""
$(TYPEDEF)
"""
mutable struct VelPose2VelPose2{T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} <: IIF.AbstractRelativeMinimize
  Zpose::Pose2Pose2{T1} #Zpose::T1
  Zvel::T2
  reuseres::Vector{Vector{Float64}}
  VelPose2VelPose2{T1,T2}() where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}()
  VelPose2VelPose2{T1,T2}(z1::T1, z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = new{T1,T2}(Pose2Pose2(z1),z2,[zeros(3) for i in 1:Threads.nthreads()])
end
VelPose2VelPose2(z1::T1, z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = VelPose2VelPose2{T1,T2}(z1, z2)

getSample(cf::CalcFactor{<:VelPose2VelPose2}, N::Int=1) = ([rand(cf.factor.Zpose.z,N);rand(cf.factor.Zvel,N)], )

function (cf::CalcFactor{<:VelPose2VelPose2})(res::AbstractVector{<:Real},
                                              meas,
                                              Xi,
                                              Xj  )
  #
  z = meas
  wxi, wxj = Xi, Xj
  # @show z, wxi, wxj
  dt = Dates.value(cf.metadata.fullvariables[2].nstime - cf.metadata.fullvariables[1].nstime)*1e-9     # roughly the intended use of userdata
  # fill!(vp2vp2.reuseres[Threads.threadid()], 0.0)
  wXjhat = SE2(wxi[1:3])*SE2(meas[1:3])
  jXjhat = SE2(wxj[1:3]) \ wXjhat
  se2vee!(cf.factor.reuseres[Threads.threadid()], jXjhat)
  # vp2vp2.Zpose(vp2vp2.reuseres[Threads.threadid()], userdata, idx, meas, Xi, Xj)

  wDXij = (wxj[4:5]-wxi[4:5])
  bDXij = TransformUtils.R(-wxi[3])*wDXij
  
  # calculate the residual
  dx = se2vee(SE2(wxi[1:3]) \ SE2(wxj[1:3]))
  res[1:3] .= cf.factor.reuseres[Threads.threadid()]
  res[4:5] .= z[4:5] .- bDXij
  res[4:5] .^= 2
  res[4:5] .+= (dx[1:2]/dt .- 0.5*(wxj[4:5] .+ wxi[4:5])).^2
  res[4:5] .= sqrt.(res[4:5])

  # res[1] = sum((cf.factor.reuseres[Threads.threadid()]).^2)
  # res[1] += sum((z[4:5] - bDXij).^2)
  # res[1] += sum((dx[1:2]/dt - 0.5*(wxj[4:5]+wxi[4:5])).^2)  # first order integration
  # res[1]

  nothing
end



function compare(a::VelPose2VelPose2, b::VelPose2VelPose2; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && RoME.compare(a.Zpose, b.Zpose)
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
  return PackedVelPose2VelPose2(convert(PackedSamplableBelief, d.Zpose.z),
                                convert(PackedSamplableBelief, d.Zvel))
end
function convert(::Type{VelPose2VelPose2}, d::PackedVelPose2VelPose2)
  posedistr = convert(SamplableBelief, d.strpose)
  veldistr = convert(SamplableBelief, d.strvel)
  return VelPose2VelPose2(posedistr, veldistr)
end
