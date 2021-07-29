
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
  z::T
end
# convenience and default constructor
Pose2Pose2() = Pose2Pose2(MvNormal([1.0; 1.0; 1.0]))
 
DFG.getManifold(::Pose2Pose2) = Manifolds.SpecialEuclidean(2)

# Towards using Manifolds.jl -- this will be easier if IIF just used vector or points, rather than current VND.val::Matrix
  # basis::DefaultOrthogonalBasis{ℝ}
  # shape::S
  # R0::R
  # p0::P
  # Pose2Pose2( z::T=MvNormal(zeros(3),LinearAlgebra.diagm([1.0;1.0;1.0])), 
  #             basis=DefaultOrthogonalBasis(),
  #             shape::S=Manifolds.ShapeSpecification(Manifolds.StaticReshaper(), Manifolds.base_manifold(Manifolds.SpecialEuclidean(2)).manifolds...),
  #             R0::R=SA[1 0; 0 1],
  #             p0::P=Manifolds.prod_point(shape, 
  #                                         (SA[0;0],R0)... 
  #                                       )
  #             ) where {T <: IIF.SamplableBelief, S, R, P} = new{T,S,R,P}( z, 
  #                                                                         basis, 
  #                                                                         shape,
  #                                                                         R0,
  #                                                                         p0  )
# end

Pose2Pose2(::UniformScaling) = Pose2Pose2() # MvNormal(zeros(3),LinearAlgebra.diagm([1.0;1.0;1.0])) )


function getSample(cf::CalcFactor{<:Pose2Pose2}, N::Int=1) 
  # Z = cf.factor.Z
  # M = getManifold(Pose2)
  # p = getPointIdentity(Pose2)
  
  # Xc = [rand(Z) for _ in 1:N]
  
  # # X = get_vector.(Ref(M), Ref(p), Xc, Ref(DefaultOrthogonalBasis()))
  # X = hat.(Ref(M), Ref(p), Xc)
  # points = exp.(Ref(M), Ref(p), X)

  # (points, )

  Xc = [rand(cf.factor.z) for _ in 1:N]
  (Xc, )
end


# function (cf::CalcFactor{<:Pose2Pose2})(meas,
#                                         wxi,
#                                         wxj  )
#   #

#     #
#     wTjhat = SE2(wxi)*SE2(meas)
#     jTjhat = SE2(wxj) \ wTjhat
#     return se2vee(jTjhat)
# end

function (cf::CalcFactor{<:Pose2Pose2})(Xc, p, q)
    M = getManifold(cf.factor)
    X = hat(M, p, Xc)

    q̂ = Manifolds.compose(M, p, exp(M, identity(M, p), X)) #for groups
    return log(M, q, q̂)
end
  
    # G = getManifold(cf.factor)
    # p0 = cf.factor.p0
    # R0 = cf.factor.R0
    # M_R = base_manifold(G)[2]
      # R0 = [1.0 0; 0 1] # cf.factor.R0
      # _t0, _w0 = [0.0; 0],  hat(M_R, R0, 0)
      # p0 = ProductRepr(zeros(2), Matrix{Float64}(I,2,2))
  
  
      # jPjhat = Manifolds.ProductRepr(jTjhat[1:2,3], jTjhat[1:2,1:2])
      # return get_coordinates(G, p0, log(G, p0, jPjhat), cf.factor.basis)
  
      # wPihat = Manifolds.prod_point(cf.factor.shape, (wxi[1:2], exp(M_R, R0, hat(M_R, R0, wxi[3])))... )
      # iPj    = Manifolds.prod_point(cf.factor.shape, (meas[1:2], exp(M_R, R0, hat(M_R, R0, meas[3])))... )
      # wPjhat = Manifolds.prod_point(cf.factor.shape, (wxj[1:2], exp(M_R, R0, hat(M_R, R0, wxj[3])))... )
  
    # wWi = get_vector(M_R, R0, wxi[3], cf.factor.basis)
    # wXi = Manifolds.ProductRepr(wxi[1:2], wWi)
  
    # iWj = get_vector(M_R, R0, meas[3], cf.factor.basis)
    # iXj = Manifolds.ProductRepr(meas[1:2], iWj)
  
    # wWj = get_vector(M_R, R0, wxj[3], cf.factor.basis)
    # wXj = Manifolds.ProductRepr(wxj[1:2], wWi)
  
    # wPihat = Manifolds.exp(G, p0, wXi)
    # iPj    = Manifolds.exp(G, p0, iXj)
    # wPjhat = Manifolds.exp(G, p0, wXj)
  
      # @show wPjpred = compose(G, wPihat, iPj)
      # @show Manifolds.affine_matrix(G, wPjpred)
  
      # @show jhatPw = Manifolds.inv(G, wPjhat)
      # @show jPjhat = compose(G, jhatPw, wPjpred)
      # @show Manifolds.affine_matrix(G, jPjhat)
  
  
    # wTi = Manifolds.affine_matrix(G, wPihat)
    # iTj = Manifolds.affine_matrix(G, iPj)
    # wTj = Manifolds.affine_matrix(G, wPjhat)
  
    # # #
    # @show jTjhat = (wTi*iTj)\wTj
    # jPjhat = Manifolds.ProductRepr(jTjhat[1:2,3], jTjhat[1:2,1:2])
    # return get_coordinates(G, p0, log(G, p0, jPjhat), cf.factor.basis)
    
      # error("now this")
  
      # @show jXj = log(G, p0, jPjhat)
      # return get_coordinates(G, p0, log(G, p0, jPjhat), cf.factor.basis)


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
  return PackedPose2Pose2(convert(PackedSamplableBelief, d.z))
end



# FIXME, rather have separate compareDensity functions
function compare(a::Pose2Pose2,b::Pose2Pose2; tol::Float64=1e-10)
  return compareDensity(a.z, b.z)
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
