# Pose3Pose3 evaluation functions





type PriorPose3 <: IncrementalInference.FunctorSingleton
    Zi::SE3
    Cov::Array{Float64,2}
    PriorPose3() = new()
    PriorPose3(st::FloatInt, sr::Float64) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]])
    PriorPose3(s::SE3, c::Array{Float64,2}) = new(s,c)
end
function getSample(p3::PriorPose3, N::Int=1)
  mv = Distributions.MvNormal(veeEuler(p3.Zi), p3.Cov)
  return (rand(mv, N),)
end
type PackedPriorPose3  <: IncrementalInference.PackedInferenceType
    vecZi::Array{Float64,1} # 0rotations, 1translation in each column
    vecCov::Array{Float64,1}
    dimc::Int64
    PackedPriorPose3() = new()
    PackedPriorPose3(x...) = new(x[1], x[2], x[3])
end
function convert(::Type{PriorPose3}, d::PackedPriorPose3)
  Zi = SE3(d.vecZi[1:3], Quaternion(d.vecZi[4],d.vecZi[5:7]))
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPose3(Zi, Cov)
end
function convert(::Type{PackedPriorPose3}, d::PriorPose3)
  # TODO -- change to
  v1 = veeQuaternion(d.Zi)
  v2 = d.Cov[:];
  return PackedPriorPose3(v1, v2, size(d.Cov,1))
end




# ------------------------------------

type PP3REUSE
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

type Pose3Pose3 <: RoME.BetweenPoses
    Zij::Distribution
    reuse::Vector{PP3REUSE}
    Pose3Pose3() = new()
    Pose3Pose3(s) = new(s, fill(PP3REUSE(), Threads.nthreads() )  )
end
function getSample(pp3::Pose3Pose3, N::Int=1)
  return (rand(pp3.Zij, N), )
end
function (pp3::Pose3Pose3)(res::Array{Float64},
      idx::Int,
      meas::Tuple,
      wXi::Array{Float64,2},
      wXj::Array{Float64,2}  )
  #
  reusethrid = pp3.reuse[Threads.threadid()]
  fastpose3pose3residual!(reusethrid, res, idx, meas, wXi, wXj)
  nothing
end


type PackedPose3Pose3 <: IncrementalInference.PackedInferenceType
  vecZij::Array{Float64,1} # 3translations, 3rotation
  vecCov::Array{Float64,1}
  dimc::Int64
  PackedPose3Pose3() = new()
  PackedPose3Pose3(x...) = new(x[1], x[2], x[3])
end
function convert(::Type{Pose3Pose3}, d::PackedPose3Pose3)
  qu = Quaternion(d.vecZij[4], d.vecZij[5:7])
  se3val = SE3(d.vecZij[1:3], qu)
  cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Pose3Pose3( MvNormal(veeEuler(se3val), cov) )
end
function convert(::Type{PackedPose3Pose3}, d::Pose3Pose3)
  val = d.Zij.μ
  se3val = SE3(val[1:3], Euler(val[4:6]...))
  v1 = veeQuaternion(se3val)
  v2 = d.Zij.Σ.mat
  return PackedPose3Pose3(v1[:], v2[:], size(v2,1) )
end




# -----------------------

type Pose3Pose3NH <: IncrementalInference.FunctorPairwiseNH
    Zij::Distribution
    nullhypothesis::Distributions.Categorical
    reuse::Vector{PP3REUSE}
    Pose3Pose3NH() = new()
    Pose3Pose3NH(s::Distribution, vh::Vector{Float64}) = new(s, Distributions.Categorical(vh), fill(PP3REUSE(), Threads.nthreads() )  )
    # Pose3Pose3NH(s::SE3, c::Array{Float64,2}, vh::Float64) = new(s,c, Distributions.Categorical([(1.0-vh);vh]),SE3(0),SE3(0),SE3(0))
    # Pose3Pose3NH(st::FloatInt, sr::Float64;vh::Float64=1.0) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]], Distributions.Categorical([(1.0-vh);vh]),SE3(0),SE3(0),SE3(0))
end
function getSample(pp3::Pose3Pose3NH, N::Int=1)
  return (rand(pp3.Zij, N), )
end
function (pp3::Pose3Pose3NH)(res::Array{Float64},
      idx::Int,
      meas::Tuple,
      wXi::Array{Float64,2},
      wXj::Array{Float64,2}  )
  #
  reusethrid = pp3.reuse[Threads.threadid()]
  fastpose3pose3residual!(reusethrid, res, idx, meas, wXi, wXj)
  nothing
end



type PackedPose3Pose3NH <: IncrementalInference.PackedInferenceType
  vecZij::Vector{Float64} # 3translations, 3rotation
  vecCov::Vector{Float64}
  dimc::Int64
  nullhypothesis::Vector{Float64}
  PackedPose3Pose3NH() = new()
  PackedPose3Pose3NH(x1::Vector{Float64},x2::Vector{Float64},x3::Int64,x4::Vector{Float64}) = new(x1, x2, x3, x4)
end

function convert(::Type{Pose3Pose3NH}, d::PackedPose3Pose3NH)
  qu = Quaternion(d.vecZij[4], d.vecZij[5:7])
  se3val = SE3(d.vecZij[1:3], qu)
  cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return Pose3Pose3NH( MvNormal(veeEuler(se3val), cov), d.nullhypothesis )
end
function convert(::Type{PackedPose3Pose3NH}, d::Pose3Pose3NH)
  val = d.Zij.μ
  se3val = SE3(val[1:3], Euler(val[4:6]...))
  v1 = veeQuaternion(se3val)
  v2 = d.Zij.Σ.mat
  return PackedPose3Pose3NH(v1[:], v2[:], size(v2,1), d.nullhypothesis.p )
end


# TODO -- stronger type safety required here
# Project all particles (columns) Xval with Z, that is for all  SE3(Xval[:,i])*Z
function projectParticles(Xval::Array{Float64,2}, Z::Distribution)
  # TODO optimize for more speed with threads and better memory management
  r,c = size(Xval)
  RES = zeros(r,c) #*cz

  ent, x = SE3(0), SE3(0)
  ENT = rand( Z, c )
  # ENT = rand( MvNormal(zeros(6), Cov), c )
  j=1
  # for j in 1:cz
  for i in 1:c
    x.R = TransformUtils.convert(SO3,Euler(Xval[4,i],Xval[5,i],Xval[6,i]))
    x.t = Xval[1:3,i]
    ent.R =  TransformUtils.convert(SO3, Euler(ENT[4:6,i]...)) # so3
    ent.t = ENT[1:3,i]
    # newval = Z*ent
    # res = x*newval
    res = x*ent
    RES[1:r,i*j] = veeEuler(res)
  end
  # end
  #
  return RES
end

⊕(Xpts::Array{Float64,2}, z::Pose3Pose3) = projectParticles(Xpts, z.Zij)
⊕(Xvert::Graphs.ExVertex, z::Pose3Pose3) = ⊕(getVal(Xvert), z)





# 3 dimensionsal evaluation functions.
# Transforms code is scattered, and will be moved to better location with stable 3D Examples
# also worth it to consolidate if other affine transforms package is available

# using Euler angles for linear sampling in product of potentials
# Gaussian model for prior
# function evalPotential(obs::PriorPose3, Xi::Array{Graphs.ExVertex,1}; N::Int64=200)
#   mu = veeEuler(obs.Zi)
#   return rand( MvNormal(mu, obs.Cov), N)
# end
#
#
# # Still limited to linear sampler, then reprojected onto ball -- TODO upgrade manifold sampler
# function evalPotential(odom::Pose3Pose3, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int64=200)
#     # rz,cz = size(odom.Zij)
#     Xval = Array{Float64,2}()
#     # implicit equation portion -- bi-directional pairwise function made explicit here
#     if Xid == Xi[1].index #odom.
#         # reverse direction
#         Z = inverse(odom.Zij)
#         Xval = getVal(Xi[2])
#     elseif Xid == Xi[2].index
#         # forward direction
#         Z = odom.Zij
#         Xval = getVal(Xi[1])
#     else
#         error("Bad evalPairwise Pose3Pose3")
#     end
#
#     return projectParticles(Xval, Z, odom.Cov)
# end
#
# # Still limited to linear sampler, then reprojected onto ball -- TODO upgrade manifold sampler
# function evalPotential(odom::Pose3Pose3NH,
#       Xi::Array{Graphs.ExVertex,1},
#       Xid::Int64;
#       N::Int64=200)
#    #
#   # rz,cz = size(odom.Zij)
#   Xval = Array{Float64,2}()
#   # implicit equation portion -- bi-directional pairwise function made explicit here
#   if Xid == Xi[1].index #odom.
#       # reverse direction
#       Z = inverse(odom.Zij)
#       Xval = getVal(Xi[2])
#       oldval = getVal(Xi[1])
#   elseif Xid == Xi[2].index
#       # forward direction
#       Z = odom.Zij
#       Xval = getVal(Xi[1])
#       oldval = getVal(Xi[2])
#   else
#       error("Bad evalPairwise Pose3Pose3")
#   end
#
#   skipdos = rand(odom.ValidHypot, N)
#   dos, donts = skipdos .== 2, skipdos .== 1
#
#   projted = Array{Float64}(6, N)
#   projted[:,dos] = projectParticles(Xval[:,dos], Z, odom.Cov)
#
#   projted[:,donts] = oldval[:,donts]
#
#   return projted
# end










#
