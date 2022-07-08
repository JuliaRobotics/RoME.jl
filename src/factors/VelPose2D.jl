

"""
$(TYPEDEF)
"""
struct VelPose2VelPose2{T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  Zpose::Pose2Pose2{T1}
  Zvel::T2
end
VelPose2VelPose2(z1::SamplableBelief, z2::SamplableBelief) = VelPose2VelPose2(Pose2Pose2(z1), z2)

getManifold(::InstanceType{VelPose2VelPose2}) = getManifold(DynPose2)

preableCache(::AbstractDFG, ::AbstractVector{Symbol}, ::VelPose2VelPose2) = zeros(3)

function getSample(cf::CalcFactor{<:VelPose2VelPose2})
    #Pose2 part
    Xc = rand(cf.factor.Zpose.Z)
    M = getManifold(Pose2)
    ϵ = getPointIdentity(Pose2)
    # ϵ = Manifolds.Identity(M)
    Xpose = hat(M, ϵ, Xc)
    #velocity part
    Xvel = rand(cf.factor.Zvel)

    return ArrayPartition(Xpose, Xvel)
end


function IIF.getMeasurementParametric(s::VelPose2VelPose2{<:MvNormal, <:MvNormal}) 

  meas = [mean(s.Zpose.z); mean(s.Zvel)]

  iΣp = invcov(s.Zpose.z)
  iΣv = invcov(s.Zvel)

  iΣ = zeros(eltype(iΣp), 5,5)

  iΣ[1:3,1:3] .= iΣp
  iΣ[4:5,4:5] .= iΣv

  return meas, iΣ
end

function (cf::CalcFactor{<:VelPose2VelPose2})(X, p, q)
  #
  pose_res = Vector{Manifolds.number_eltype(X)}(undef, 3)
  #Pose2 part
  M1 = getManifold(Pose2)
  X1 = submanifold_component(X,1) 
  p1 = submanifold_component(p,1) 
  q1 = submanifold_component(q,1) 
  ϵ1 = identity_element(M1, p1)
  q̂1 = Manifolds.compose(M1, p1, exp(M1, ϵ1, X1))
  vee!(M1, pose_res, q1, log(M1, q1, q̂1))
  
  #velocity part
  dt = Dates.value(cf.metadata.fullvariables[2].nstime - cf.metadata.fullvariables[1].nstime)*1e-9
  X2 = submanifold_component(X,2)
  p2 = submanifold_component(p,2)
  q2 = submanifold_component(q,2)
  # bDXij = TransformUtils.R(-wxi[3])*wDXij
  bDXij = transpose(submanifold_component(p1,2))*(q2 .- p2)

  Xpq = log(M1, ϵ1, Manifolds.compose(M1, Manifolds.inv(M1, p1), q1))
  dx = Vector{Manifolds.number_eltype(X)}(undef, 3)
  vee!(M1, dx, ϵ1, Xpq)
  # calculate the residual
  res_vel = (X2 .- bDXij).^2 .+ (view(dx, 1:2)/dt .- 0.5*(p2 .+ q2)).^2
  res_vel = sqrt.(res_vel)

  return [pose_res; res_vel]
end
#=
function (cf::CalcFactor{<:VelPose2VelPose2})(meas,
                                              _Xi,
                                              _Xj  )
  #
  M = getManifold(DynPose2)
  Xi = vee(M, _Xi, log(M, identity_element(M, _Xi), _Xi))
  Xj = vee(M, _Xj, log(M, identity_element(M, _Xj), _Xj))

  #FIXME JT - Createing new res for simplicity, it may not hold up well though
  res = Vector{eltype(Xi)}(undef, 5)
  
  z = meas
  wxi, wxj = Xi, Xj
  # @show z, wxi, wxj
  dt = Dates.value(cf.metadata.fullvariables[2].nstime - cf.metadata.fullvariables[1].nstime)*1e-9     # roughly the intended use of userdata
  # fill!(vp2vp2.reuseres[Threads.threadid()], 0.0)
  wXjhat = SE2(wxi[1:3])*SE2(meas[1:3])
  jXjhat = SE2(wxj[1:3]) \ wXjhat

  #FIXME this does not work with parametric
  # se2vee!(cf.factor.reuseres[Threads.threadid()], jXjhat)
  #FIXME cf.factor.reuseres has type issues with parametric
  se2vee!(res, jXjhat)

  # vp2vp2.Zpose(vp2vp2.reuseres[Threads.threadid()], userdata, idx, meas, Xi, Xj)

  wDXij = (wxj[4:5]-wxi[4:5])
  bDXij = TransformUtils.R(-wxi[3])*wDXij
  
  # calculate the residual
  dx = se2vee(SE2(wxi[1:3]) \ SE2(wxj[1:3]))
  #FIXME cf.factor.reuseres has type issues with parametric
  # res[1:3] .= cf.factor.reuseres[Threads.threadid()]
  res[4:5] .= z[4:5] .- bDXij
  res[4:5] .^= 2
  res[4:5] .+= (dx[1:2]/dt .- 0.5*(wxj[4:5] .+ wxi[4:5])).^2
  res[4:5] .= sqrt.(res[4:5])

  # res[1] = sum((cf.factor.reuseres[Threads.threadid()]).^2)
  # res[1] += sum((z[4:5] - bDXij).^2)
  # res[1] += sum((dx[1:2]/dt - 0.5*(wxj[4:5]+wxi[4:5])).^2)  # first order integration
  # res[1]

  return res
end
=#


function compare(a::VelPose2VelPose2, b::VelPose2VelPose2; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && RoME.compare(a.Zpose, b.Zpose)
  TP = TP && RoME.compareDensity(a.Zvel, b.Zvel)
  # TP = TP && norm(a.reuseres - b.reuseres) < tol
  return TP
end


"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedVelPose2VelPose2 <: AbstractPackedFactor
  Zpose::PackedSamplableBelief
  Zvel::PackedSamplableBelief
end

function convert(::Type{PackedVelPose2VelPose2}, d::VelPose2VelPose2)
  return PackedVelPose2VelPose2(convert(PackedSamplableBelief, d.Zpose.Z),
                                convert(PackedSamplableBelief, d.Zvel))
end
function convert(::Type{VelPose2VelPose2}, d::PackedVelPose2VelPose2)
  posediZ = convert(SamplableBelief, d.Zpose)
  veldiZ = convert(SamplableBelief, d.Zvel)
  return VelPose2VelPose2(posediZ, veldiZ)
end
