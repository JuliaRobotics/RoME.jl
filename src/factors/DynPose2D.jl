# 2D SLAM with velocity states


"""
$(TYPEDEF)
"""
mutable struct DynPose2VelocityPrior{T1,T2} <: IncrementalInference.AbstractPrior where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief}
  Zpose::T1
  Zvel::T2
end
DynPose2VelocityPrior(z1::T1,z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = DynPose2VelocityPrior{T1,T2}(z1,z2)

DFG.getManifold(::DynPose2VelocityPrior) = getManifold(DynPose2)

function getSample(cf::CalcFactor{<:DynPose2VelocityPrior})
  Zpose = cf.factor.Zpose
  Zvel = cf.factor.Zvel
  p = getPointIdentity(DynPose2())
  M = cf.manifold # getManifold(cf.factor)
  
  Xc = [rand(Zpose);rand(Zvel)]
  
  # X = get_vector.(Ref(M), Ref(p), Xc, Ref(DefaultOrthogonalBasis()))
  X = hat(M, p, Xc)
  points = exp(M, p, X)

  return points
end

function IIF.getMeasurementParametric(s::DynPose2VelocityPrior{<:MvNormal, <:MvNormal}) 

  meas = [mean(s.Zpose); mean(s.Zvel)]

  iΣp = invcov(s.Zpose)
  iΣv = invcov(s.Zvel)

  iΣ = zeros(eltype(iΣp), 5,5)

  iΣ[1:3,1:3] .= iΣp
  iΣ[4:5,4:5] .= iΣv

  return meas, iΣ
end

function (cf::CalcFactor{<:DynPose2VelocityPrior})(meas, X)
  # FIXME update Manifolds.jl usage
  #pose part, reused from PriorPose2
  iXihat = SE2(meas[1:3]) \ SE2(X[1:3])	
  res_pose = se2vee(iXihat)	

  #velocity part, normal prior
  res_vel = meas[4:5] .- X[4:5] 

  return [res_pose; res_vel]
end

"""
$(TYPEDEF)
"""
mutable struct DynPose2Pose2{T <: IIF.SamplableBelief} <: IIF.AbstractRelativeRoots
  Zpose::Pose2Pose2{T} #Zpose::T1
  partial::Tuple{Int,Int,Int}
end

preableCache(::AbstractDFG, ::AbstractVector{Symbol}, ::DynPose2Pose2) = zeros(5)

getSample(cf::CalcFactor{<:DynPose2Pose2}) = rand(cf.factor.Zpose.z)

function (cf::CalcFactor{<:DynPose2Pose2})( meas,
                                            wXi,
                                            wXj  )
  #
  # FIXME update to Manifolds.jl usage
  # cf.factor.Zpose(res, meas, wXi, wXj)
  wXjhat = SE2(wXi)*SE2(meas)
  jXjhat = SE2(wXj) \ wXjhat
  return se2vee(jXjhat)
  # z = meas[1][:,idx]
  # wxi, wxj = wXi[:,idx], wXj[:,idx]
  # dt = (userdata.variableuserdata[2].ut - userdata.variableuserdata[1].ut)*1e-6
  # wpj = ( wxi[1:2]+dt*wxi[4:5] + z[1:2] )
  # thetaj = se2vee(SE2([0;0;wxi[3]])*SE2([0;0;z[3]]))[3]
  # res[1:3] = se2vee( SE2(wxj[1:3])\SE2([wpj;thetaj]) )
  # # res[1:2] = (wXj[1:2] - (wXi[1:2]+dt*wXi[4:5])+z[1:2])
  # res[4:5] = z[4:5] - (wxj[4:5] -  wxi[4:5])
  # nothing
end




function compare(a::DynPose2VelocityPrior, b::DynPose2VelocityPrior)::Bool
  RoME.compareDensity(a.Zpose, b.Zpose) && RoME.compareDensity(a.Zvel, b.Zvel)
end




## Packing types


"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedDynPose2VelocityPrior <: AbstractPackedFactor
  strpose::PackedSamplableBelief
  strvel::PackedSamplableBelief
end

function convert(::Type{PackedDynPose2VelocityPrior}, d::DynPose2VelocityPrior)
  return PackedDynPose2VelocityPrior(convert(PackedSamplableBelief, d.Zpose),convert(PackedSamplableBelief, d.Zvel))
end
function convert(::Type{DynPose2VelocityPrior}, d::PackedDynPose2VelocityPrior)
  posedistr = convert(SamplableBelief, d.strpose)
  veldistr = convert(SamplableBelief, d.strvel)
  return DynPose2VelocityPrior(posedistr, veldistr)
end


"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedDynPose2Pose2 <: AbstractPackedFactor
  Z::PackedSamplableBelief
end

function convert(::Type{PackedDynPose2Pose2}, d::DynPose2Pose2)
  return PackedDynPose2Pose2(convert(PackedSamplableBelief, d.Zpose.Z))
end
function convert(::Type{DynPose2Pose2}, d::PackedDynPose2Pose2)
  return DynPose2Pose2(convert(SamplableBelief, d.Z))
end




## ===========================================================================================



"""
$(TYPEDEF)
"""
Base.@kwdef struct DynPose2DynPose2{T <: IIF.SamplableBelief} <: AbstractRelativeRoots
  Z::T = MvNormal(zeros(5), diagm([0.01;0.01;0.001;0.1;0.1].^2))
end
preambleCache(::AbstractDFG, ::AbstractVector{<:DFGVariable}, ::DynPose2DynPose2) = zeros(5)

# FIXME ON FIRE, must update to new Manifolds style factors
getManifold(::DynPose2DynPose2) = SE2E2_Manifold # not fully impl manifold yet
# FIXME, should produce tangents, not coordinates.
getSample(cf::CalcFactor{<:DynPose2DynPose2}) = rand(cf.factor.Z)
function (cf::CalcFactor{<:DynPose2DynPose2})(meas,
                                              wXi,
                                              wXj  )
  #
  # FIXME update to Manifolds.jl usage

  # vp2vp2.Zpose(res, userdata, idx, meas, Xi, Xj)
    # wXjhat = SE2(wxi[1:3,idx])*SE2(meas[1][1:3,idx])
    # jXjhat = SE2(wxj[1:3,idx]) \ wXjhat
    # se2vee!(res, jXjhat)
  z = meas
  wxi, wxj = wXi, wXj
  dt = Dates.value(cf.fullvariables[2].nstime - cf.fullvariables[1].nstime)*1e-9  
  wpj = ( wxi[1:2]+dt*wxi[4:5] + z[1:2] )
  thetaj = se2vee(SE2([0;0;wxi[3]])*SE2([0;0;z[3]]))[3]
  res13 = se2vee( SE2(wxj[1:3])\SE2([wpj;thetaj]) )
  # res[1:2] = (wXj[1:2] - (wXi[1:2]+dt*wXi[4:5])+z[1:2])
  res45 = z[4:5] - (wxj[4:5] -  wxi[4:5])
  return [res13;res45]
end



function compare(a::DynPose2DynPose2, b::DynPose2DynPose2; tol::Float64=1e-10)::Bool
  TP = true
  TP = TP && RoME.compareDensity(a.Z, b.Z)
  # TP = TP && norm(a.reuseres - b.reuseres) < tol
  return TP
end


"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedDynPose2DynPose2 <: AbstractPackedFactor
  Z::PackedSamplableBelief
end

function convert(::Type{PackedDynPose2DynPose2}, d::DynPose2DynPose2)
  return PackedDynPose2DynPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{DynPose2DynPose2}, d::PackedDynPose2DynPose2)
  return DynPose2DynPose2(convert(SamplableBelief, d.Z))
end





#
