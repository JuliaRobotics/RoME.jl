# 2D SLAM with velocity states

export   DynPose2, DynPose2VelocityPrior, PackedDynPose2VelocityPrior, DynPose2Pose2, PackedDynPose2Pose2



## ==================================================================================================
## Hack to be removed or updated
## FIXME ME ON FIRE use the ManifoldsBase.jl prescribed interface method instead:
##   https://juliamanifolds.github.io/Manifolds.jl/stable/examples/manifold.html#manifold-tutorial
## ==================================================================================================


import ApproxManifoldProducts: coords, uncoords, getPointsManifold, _makeVectorManifold

export SE2E2_Manifold


# this is a hack and not fully implemented as per: 
struct _SE2E2 <: ManifoldsBase.Manifold{ManifoldsBase.ℝ} end

const SE2E2_Manifold = _SE2E2()

AMP.coords(::Type{<:typeof(SE2E2_Manifold)}, p::ProductRepr) = [p.parts[1][1], p.parts[1][2], atan(p.parts[2][2,1],p.parts[2][1,1]), p.parts[3][1], p.parts[3][2]]

function AMP.uncoords(::typeof(SE2E2_Manifold), p::AbstractVector{<:Real})
  α = p[3]
  ProductRepr(([p[1], p[2]]), [cos(α) -sin(α); sin(α) cos(α)], ([p[4], p[5]]))
end

function AMP.getPointsManifold(mkd::ManifoldKernelDensity{M}) where {M <: typeof(SE2E2_Manifold)}
  data_ = getPoints(mkd.belief)
  [uncoords(mkd.manifold, view(data_, :, i)) for i in 1:size(data_,2)]
end

function Statistics.mean(::typeof(SE2E2_Manifold), pts::AbstractVector)
  se2_ = (d->ProductRepr(d.parts[1], d.parts[2])).(pts)
  mse2 = mean(SpecialEuclidean(2), se2_)
  e2_ = (d->ProductRepr(d.parts[3])).(pts)
  me2 = mean(Euclidean(2), e2_)
  ProductRepr(mse2.parts[1], mse2.parts[2], me2.parts[1])
end

AMP._makeVectorManifold(::M, prr::ProductRepr) where {M <: typeof(SE2E2_Manifold)} = coords(M, prr)



## ==================================================================================================



"""
$(TYPEDEF)
"""
mutable struct DynPose2VelocityPrior{T1,T2} <: IncrementalInference.AbstractPrior where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief}
  Zpose::T1
  Zvel::T2
  DynPose2VelocityPrior{T1,T2}() where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = new{T1,T2}()
  DynPose2VelocityPrior{T1,T2}(z1::T1,z2::T2) where {T1 <: IIF.SamplableBelief,T2 <: IIF.SamplableBelief} = new{T1,T2}(z1,z2)
end
DynPose2VelocityPrior(z1::T1,z2::T2) where {T1 <: IIF.SamplableBelief, T2 <: IIF.SamplableBelief} = DynPose2VelocityPrior{T1,T2}(z1,z2)

getSample(cf::CalcFactor{<:DynPose2VelocityPrior}, N::Int=1) = ([rand(cf.factor.Zpose,N);rand(cf.factor.Zvel,N)], )


function IIF.getParametricMeasurement(s::DynPose2VelocityPrior{<:MvNormal, <:MvNormal}) 

  meas = [mean(s.Zpose); mean(s.Zvel)]

  iΣp = invcov(s.Zpose)
  iΣv = invcov(s.Zvel)

  iΣ = zeros(eltype(iΣp), 5,5)

  iΣ[1:3,1:3] .= iΣp
  iΣ[4:5,4:5] .= iΣv

  return meas, iΣ
end

function (cf::CalcFactor{<:DynPose2VelocityPrior})(meas, X)
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
  # reuseres::Vector{Float64}
  partial::Tuple{Int,Int,Int}
  DynPose2Pose2{T}() where {T <: IIF.SamplableBelief} = new{T}()
  DynPose2Pose2{T}(z1::T) where {T <: IIF.SamplableBelief} = new{T}(Pose2Pose2(z1), (1,2,3))
end
DynPose2Pose2(z1::T) where {T <: IIF.SamplableBelief} = DynPose2Pose2{T}(z1)

getSample(cf::CalcFactor{<:DynPose2Pose2}, N::Int=1) = (rand(cf.factor.Zpose.z,N), )

function (cf::CalcFactor{<:DynPose2Pose2})( meas,
                                            wXi,
                                            wXj  )
  #
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
mutable struct PackedDynPose2VelocityPrior <: IncrementalInference.PackedInferenceType
  strpose::AbstractString
  strvel::AbstractString
  PackedDynPose2VelocityPrior() = new()
  PackedDynPose2VelocityPrior(z1::AS, z2::AS) where {AS <: AbstractString} = new(z1, z2)
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




## ===========================================================================================



"""
$(TYPEDEF)
"""
mutable struct DynPose2DynPose2{T <: IIF.SamplableBelief} <: AbstractRelativeRoots
  Z::T
  reuseres::Vector{Vector{Float64}}
  DynPose2DynPose2{T}() where {T <: IIF.SamplableBelief} = new{T}()
  DynPose2DynPose2{T}(z1::T) where {T <: IIF.SamplableBelief} = new{T}(z1,[zeros(5) for i in 1:Threads.nthreads()])
end
DynPose2DynPose2(z1::T=MvNormal(zeros(5), diagm([0.01;0.01;0.001;0.1;0.1].^2))) where {T <: IIF.SamplableBelief} = DynPose2DynPose2{T}(z1)

getManifolds(::Type{DynPose2DynPose2}) = (:Euclid,:Euclid,:Circular,:Euclid,:Euclid)
getManifolds(::DynPose2DynPose2) = getManifolds(DynPose2DynPose2)


getSample(cf::CalcFactor{<:DynPose2DynPose2}, N::Int=1) = (rand(cf.factor.Z, N), )

function (cf::CalcFactor{<:DynPose2DynPose2})(meas,
                                              wXi,
                                              wXj  )
  #
  # vp2vp2.Zpose(res, userdata, idx, meas, Xi, Xj)
    # wXjhat = SE2(wxi[1:3,idx])*SE2(meas[1][1:3,idx])
    # jXjhat = SE2(wxj[1:3,idx]) \ wXjhat
    # se2vee!(res, jXjhat)
  z = meas
  wxi, wxj = wXi, wXj
  dt = Dates.value(cf.metadata.fullvariables[2].nstime - cf.metadata.fullvariables[1].nstime)*1e-9  
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
  TP = TP && norm(a.reuseres - b.reuseres) < tol
  return TP
end


"""
$(TYPEDEF)
"""
mutable struct PackedDynPose2DynPose2 <: IncrementalInference.PackedInferenceType
  Z::String
  PackedDynPose2DynPose2() = new()
  PackedDynPose2DynPose2(z1::AbstractString) = new(z1)
end

function convert(::Type{PackedDynPose2DynPose2}, d::DynPose2DynPose2)
  return PackedDynPose2DynPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{DynPose2DynPose2}, d::PackedDynPose2DynPose2)
  return DynPose2DynPose2(convert(SamplableBelief, d.Z))
end





#
