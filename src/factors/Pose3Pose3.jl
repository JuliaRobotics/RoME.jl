# Pose3Pose3 evaluation functions

export Pose3Pose3, PackedPose3Pose3

# ------------------------------------

"""
$(TYPEDEF)
"""
mutable struct PP3REUSE
  wTi::SE3
  wTj::SE3
  iTi::SE3
  PP3REUSE() = new(SE3(0),SE3(0),SE3(0))
end

function fastpose3pose3residual!(reusethrid::PP3REUSE,
                                 res::Array{Float64},
                                 idx::Int,
                                 meas::Tuple,
                                 wXi::Array{Float64,2},
                                 wXj::Array{Float64,2}  )
  #
  reusethrid.wTi.t[1:3] = wXi[1:3,idx]
  TransformUtils.convert!(reusethrid.wTi.R, Euler(wXi[4,idx],wXi[5,idx],wXi[6,idx]))
  reusethrid.wTj.t[1:3] = wXj[1:3,idx]
  TransformUtils.convert!(reusethrid.wTj.R, Euler(wXj[4,idx],wXj[5,idx],wXj[6,idx]))

  # TODO -- convert to in place convert! functions, many speed-ups possible here
  jTi = SE3( matrix(reusethrid.wTj)\matrix(reusethrid.wTi) )
  # also wasted memory here, should operate directly on iTi and not be assigning new memory
  reusethrid.iTi = (SE3(meas[1][1:3,idx],Euler(meas[1][4:6,idx]...)) * jTi)
  res[:] = veeEuler(reusethrid.iTi)
  nothing
end

"""
$(TYPEDEF)

Rigid transform factor between two Pose3 compliant variables.
"""
mutable struct Pose3Pose3{T <: IIF.SamplableBelief} <: AbstractRelativeFactor
    Zij::T
    reuse::Vector{PP3REUSE}
    Pose3Pose3{T}() where T = new{T}()
    Pose3Pose3{T}(s::T) where {T <: SamplableBelief} = new{T}(s, PP3REUSE[PP3REUSE() for i in 1:Threads.nthreads()]  )
end
Pose3Pose3(z::T=MvNormal(zeros(6),LinearAlgebra.diagm([0.01*ones(3);0.0001*ones(3)]))) where {T <: IIF.SamplableBelief} = Pose3Pose3{T}(z)

function getSample(pp3::Pose3Pose3, N::Int=1)
  return (rand(pp3.Zij, N), )
end
function (pp3::Pose3Pose3)(res::Array{Float64},
                           userdata::FactorMetadata,
                           idx::Int,
                           meas::Tuple,
                           wXi::Array{Float64,2},
                           wXj::Array{Float64,2}  )
  #
  reusethrid = pp3.reuse[Threads.threadid()]
  fastpose3pose3residual!(reusethrid, res, idx, meas, wXi, wXj)
  nothing
end

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
  # qu = Quaternion(d.vecZij[4], d.vecZij[5:7])
  # se3val = SE3(d.vecZij[1:3], qu)
  # cov = reshapeVec2Mat(d.vecCov, d.dimc)
  # return Pose3Pose3( MvNormal(veeEuler(se3val), cov) )
  return Pose3Pose3( extractdistribution(packed.Zij) )
end
function convert(::Type{PackedPose3Pose3}, obj::Pose3Pose3)
  # val = d.Zij.μ
  # se3val = SE3(val[1:3], Euler(val[4:6]...))
  # v1 = veeQuaternion(se3val)
  # v2 = d.Zij.Σ.mat
  # return PackedPose3Pose3(v1[:], v2[:], size(v2,1) )
  return PackedPose3Pose3( string(obj.Zij) )
end









#
