# Pose3Pose3 evaluation functions

export Pose3Pose3, PackedPose3Pose3

# ------------------------------------
"""
$(TYPEDEF)

Rigid transform factor between two Pose3 compliant variables.
"""
struct Pose3Pose3{T <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  z::T
end
# convenience and default constructor
Pose3Pose3() = Pose3Pose3(MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)])))

DFG.getManifold(::Pose3Pose3) = Manifolds.SpecialEuclidean(3)

Pose3Pose3(::UniformScaling) = Pose3Pose3() # MvNormal(zeros(3),LinearAlgebra.diagm([1.0;1.0;1.0])) )

function getSample(cf::CalcFactor{<:Pose3Pose3}, N::Int=1) 
  
  Xc = rand(cf.factor.z)
  #NOTE be carefull using this as template as it will not work in general for all manifolds
  M = getManifold(Pose3)
  ϵ = getPointIdentity(Pose3)
  X = hat(M, ϵ, Xc)
  # return a vector
  return (X, )
end

function (cf::CalcFactor{<:Pose3Pose3})(X, p, q)
    M = getManifold(Pose3)
    q̂ = Manifolds.compose(M, p, exp(M, identity_element(M, p), X)) #for groups
    #TODO allocalte for vee! see Manifolds #412, fix for AD
    Xc = zeros(6)
    vee!(M, Xc, q, log(M, q, q̂))
    return Xc
end

# TODO Deprecate

# """
# $(TYPEDEF)
# """
# mutable struct PP3REUSE
#   wTi::SE3
#   wTj::SE3
#   iTi::SE3
#   PP3REUSE() = new(SE3(0),SE3(0),SE3(0))
# end

# function fastpose3pose3residual!( reusethrid::PP3REUSE,
#                                   meas,
#                                   wXi,
#                                   wXj  )
#   #
#   reusethrid.wTi.t[1:3] = wXi[1:3]
#   TransformUtils.convert!(reusethrid.wTi.R, Euler(wXi[4],wXi[5],wXi[6]))
#   reusethrid.wTj.t[1:3] = wXj[1:3]
#   TransformUtils.convert!(reusethrid.wTj.R, Euler(wXj[4],wXj[5],wXj[6]))

#   # TODO -- convert to in place convert! functions, many speed-ups possible here
#   jTi = SE3( matrix(reusethrid.wTj)\matrix(reusethrid.wTi) )
#   # also wasted memory here, should operate directly on iTi and not be assigning new memory
#   reusethrid.iTi = (SE3(meas[1:3],Euler(meas[4:6]...)) * jTi)
#   return veeEuler(reusethrid.iTi)
# end

# """
# $(TYPEDEF)

# Rigid transform factor between two Pose3 compliant variables.
# """
# mutable struct Pose3Pose3{T <: IIF.SamplableBelief} <: AbstractRelativeRoots
#     Zij::T
#     reuse::Vector{PP3REUSE}
#     Pose3Pose3{T}() where T = new{T}()
#     Pose3Pose3{T}(s::T) where {T <: SamplableBelief} = new{T}(s, PP3REUSE[PP3REUSE() for i in 1:Threads.nthreads()]  )
# end
# Pose3Pose3(z::T=MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)]))) where {T <: IIF.SamplableBelief} = Pose3Pose3{T}(z)

# function getSample(cf::CalcFactor{<:Pose3Pose3}, N::Int=1)
#   return (rand(cf.factor.Zij, N), )
# end
# function (cf::CalcFactor{<:Pose3Pose3})(meas,
#                                         wXi,
#                                         wXj  )
#   #
#   reusethrid = cf.factor.reuse[Threads.threadid()]
#   return fastpose3pose3residual!(reusethrid, meas, wXi, wXj)
# end

#TODO Serialization

"""
$(TYPEDEF)

Serialization type for `Pose3Pose3`.
"""
mutable struct PackedPose3Pose3 <: IncrementalInference.PackedInferenceType
  Zij::String
  PackedPose3Pose3() = new()
  PackedPose3Pose3(x::AbstractString) = new(x)
end
function convert(::Type{Pose3Pose3}, packed::PackedPose3Pose3)
  return Pose3Pose3( convert(SamplableBelief, packed.Zij) )
end
function convert(::Type{PackedPose3Pose3}, obj::Pose3Pose3)
  return PackedPose3Pose3( convert(PackedSamplableBelief, obj.z) )
end









#
