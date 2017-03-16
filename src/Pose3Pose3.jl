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

# no longer needed -- processed by multiple dispatch
# function convert(::Type{PackedFunctionNodeData{PackedPriorPose3}}, d::FunctionNodeData{PriorPose3})
#   return PackedFunctionNodeData{PackedPriorPose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PackedPriorPose3, d.fnc))
# end
# function convert(::Type{FunctionNodeData{PriorPose3}}, d::PackedFunctionNodeData{PackedPriorPose3})
#   return FunctionNodeData{PriorPose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PriorPose3, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{PriorPose3})
#   return convert(PackedFunctionNodeData{PackedPriorPose3}, d)
# end
# function FNDdecode(d::PackedFunctionNodeData{PackedPriorPose3})
#   return convert(FunctionNodeData{PriorPose3}, d)
# end



# ------------------------------------


type Pose3Pose3 <: RoME.BetweenPoses # IncrementalInference.FunctorPairwise
    Zij::SE3 # 3translations, 3exponential param rotation, iZj
    Cov::Array{Float64,2}
    reuseTent::SE3
    reusewTi::SE3
    reusewTj::SE3
    reuseiTi::SE3
    Pose3Pose3() = new()
    Pose3Pose3(s::SE3, c::Array{Float64,2}) = new(s,c,SE3(0),SE3(0),SE3(0),SE3(0))
    Pose3Pose3(st::FloatInt, sr::Float64) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]], SE3(0),SE3(0),SE3(0),SE3(0))
end
function (pp3::Pose3Pose3)(res::Array{Float64}, idx::Int, meas::Tuple{Array{Float64,2}}, wXi::Array{Float64,2}, wXj::Array{Float64,2})
  pp3.reusewTi.t = wXi[1:3,idx]
  TransformUtils.convert!(pp3.reusewTi.R, Euler(wXi[4,idx],wXi[5,idx],wXi[6,idx]))
  pp3.reusewTj.t = wXj[1:3,idx]
  TransformUtils.convert!(pp3.reusewTj.R, Euler(wXj[4,idx],wXj[5,idx],wXj[6,idx]))

  # TODO -- convert to in place convert! functions, many speed-ups possible here
  # jTi = A_invB(SE3(0), pp3.reusewTj)*pp3.reusewTi
  # jTi = inverse(pp3.reusewTj)*pp3.reusewTi

  pp3.reuseTent.R, pp3.reuseTent.t = TransformUtils.convert(SO3, Euler(meas[1][4:6,idx]...)), meas[1][1:3,idx]
  jTi = SE3( matrix(pp3.reusewTj)\matrix(pp3.reusewTi) )
  # pp3.reuseTent.R, pp3.reuseTent.t = TransformUtils.convert(SO3, so3(meas[1][4:6,idx])), meas[1][1:3,idx]
  # pp3.reuseiTi = (pp3.Zij*pp3.reuseTent) * jTi
  pp3.reuseiTi = (pp3.reuseTent) * jTi
  res[:] = veeEuler(pp3.reuseiTi)
  nothing
end
function getSample(pp3::Pose3Pose3, N::Int=1)
  # this could be much better if we can operate with array of manifolds instead
  mv = Distributions.MvNormal(veeEuler(pp3.Zij), pp3.Cov)
  return (rand(mv, N),)
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
  return Pose3Pose3(SE3(d.vecZij[1:3], qu),
                    reshapeVec2Mat(d.vecCov, d.dimc))
end
function convert(::Type{PackedPose3Pose3}, d::Pose3Pose3)
  v1 = veeQuaternion(d.Zij)
  v2 = d.Cov[:];
  return PackedPose3Pose3(v1,v2,size(d.Cov,1) )
end


# function convert(::Type{PackedFunctionNodeData{PackedPose3Pose3}}, d::FunctionNodeData{Pose3Pose3})
#   return PackedFunctionNodeData{PackedPose3Pose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PackedPose3Pose3, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Pose3Pose3}}, d::PackedFunctionNodeData{PackedPose3Pose3})
#   return FunctionNodeData{Pose3Pose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(Pose3Pose3, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{Pose3Pose3})
#   return convert(PackedFunctionNodeData{PackedPose3Pose3}, d)
# end
# function FNDdecode(d::PackedFunctionNodeData{PackedPose3Pose3})
#   return convert(FunctionNodeData{Pose3Pose3}, d)
# end



# -----------------------



type Pose3Pose3NH <: IncrementalInference.FunctorPairwise
    Zij::SE3 # 3translations, 3exponential param rotation
    Cov::Array{Float64,2}
    ValidHypot::Distributions.Categorical
    reusewTi::SE3
    reusewTj::SE3
    reuseiTi::SE3
    Pose3Pose3NH() = new()
    Pose3Pose3NH(s::SE3, c::Array{Float64,2}, vh::Vector{Float64}) = new(s,c, Distributions.Categorical(vh),SE3(0),SE3(0),SE3(0)  )
    Pose3Pose3NH(s::SE3, c::Array{Float64,2}, vh::Float64) = new(s,c, Distributions.Categorical([(1.0-vh);vh]),SE3(0),SE3(0),SE3(0))
    Pose3Pose3NH(st::FloatInt, sr::Float64;vh::Float64=1.0) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]], Distributions.Categorical([(1.0-vh);vh]),SE3(0),SE3(0),SE3(0))
end
function getSample(pp3::Pose3Pose3NH, N::Int=1)
  # this could be much better if we can operate with array of manifolds instead
  mv = Distributions.MvNormal(veeEuler(pp3.Zij), pp3.Cov)
  return (rand(mv, N), rand(pp3.ValidHypot, N)-1)
end
function (pp3::Pose3Pose3NH)(res::Array{Float64},
      idx::Int,
      meas::Tuple,
      wXi::Array{Float64,2},
      wXj::Array{Float64,2}  )
  #
  if meas[2][idx] == 0
    # null hypothesis
    res[1:6] = 1e-10*randn(6)
  else
    pp3.reusewTi.t = wXi[1:3,idx]
    TransformUtils.convert!(pp3.reusewTi.R, Euler(wXi[4,idx],wXi[5,idx],wXi[6,idx]))
    pp3.reusewTj.t = wXj[1:3,idx]
    TransformUtils.convert!(pp3.reusewTj.R, Euler(wXj[4,idx],wXj[5,idx],wXj[6,idx]))

    # TODO -- convert to in place convert! functions, many speed-ups possible here
    jTi = SE3( matrix(pp3.reusewTj)\matrix(pp3.reusewTi) )
    pp3.reuseiTi = (SE3(meas[1][1:3,idx],Euler(meas[1][4:6,idx]...)) * jTi)
    res[:] = veeEuler(pp3.reuseiTi)
  end
  nothing
end

# using NLsolve
# function f(x, res)
#   res[1] = 1e-10*randn()
#   nothing
# end
# r = nlsolve(f, [1.0])
# @show r


type PackedPose3Pose3NH <: IncrementalInference.PackedInferenceType
  vecZij::Array{Float64,1} # 3translations, 3rotation
  vecCov::Array{Float64,1}
  dimc::Int64
  ValidHypot::Vector{Float64}
  PackedPose3Pose3NH() = new()
  PackedPose3Pose3NH(x...) = new(x[1], x[2], x[3], x[4])
end
# function convert(::Type{Pose3Pose3NH}, d::PackedPose3Pose3NH)
#   Eu = Euler(d.vecZij[4:6]...)
#   return Pose3Pose3NH(SE3(d.vecZij[1:3], Eu),
#                     reshapeVec2Mat(d.vecCov, d.dimc), d.ValidHypot)
# end
# function convert(::Type{PackedPose3Pose3NH}, d::Pose3Pose3NH)
#   v1 = veeEuler(d.Zij)
#   v2 = d.Cov[:];
#   return PackedPose3Pose3NH(v1,v2,size(d.Cov,1), d.ValidHypot.p )
# end
function convert(::Type{Pose3Pose3NH}, d::PackedPose3Pose3NH)
  qu = Quaternion(d.vecZij[4], d.vecZij[5:7])
  return Pose3Pose3NH(SE3(d.vecZij[1:3], qu),
                    reshapeVec2Mat(d.vecCov, d.dimc),
                    d.ValidHypot )
end
function convert(::Type{PackedPose3Pose3NH}, d::Pose3Pose3NH)
  v1 = veeQuaternion(d.Zij)
  v2 = d.Cov[:]
  return PackedPose3Pose3NH(v1,v2,size(d.Cov,1), d.ValidHypot.p )
end


# Project all particles (columns) Xval with Z, that is for all  SE3(Xval[:,i])*Z
function projectParticles(Xval::Array{Float64,2}, Z::SE3, Cov::Array{Float64,2})
  # TODO optimize for more speed with threads and better memory management
  r,c = size(Xval)
  RES = zeros(r,c) #*cz

  ent, x = SE3(0), SE3(0)
  ENT = rand( MvNormal(zeros(6), Cov), c )
  j=1
  # for j in 1:cz
    for i in 1:c
      x.R = TransformUtils.convert(SO3,Euler(Xval[4,i],Xval[5,i],Xval[6,i]))
      x.t = Xval[1:3,i]
      ent.R =  TransformUtils.convert(SO3, so3(ENT[4:6,i]))
      ent.t = ENT[1:3,i]
      newval = Z*ent
      res = x*newval
      RES[1:r,i*j] = veeEuler(res)
    end
  # end
  #
  return RES
end

⊕(Xpts::Array{Float64,2}, z::Pose3Pose3) = projectParticles(Xpts, z.Zij, z.Cov)
⊕(Xvert::Graphs.ExVertex, z::Pose3Pose3) = ⊕(getVal(Xvert), z)





# 3 dimensionsal evaluation functions.
# Transforms code is scattered, and will be moved to better location with stable 3D Examples
# also worth it to consolidate if other affine transforms package is available

# using Euler angles for linear sampling in product of potentials
# Gaussian model for prior
function evalPotential(obs::PriorPose3, Xi::Array{Graphs.ExVertex,1}; N::Int64=200)
  mu = veeEuler(obs.Zi)
  return rand( MvNormal(mu, obs.Cov), N)
end


# Still limited to linear sampler, then reprojected onto ball -- TODO upgrade manifold sampler
function evalPotential(odom::Pose3Pose3, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int64=200)
    # rz,cz = size(odom.Zij)
    Xval = Array{Float64,2}()
    # implicit equation portion -- bi-directional pairwise function made explicit here
    if Xid == Xi[1].index #odom.
        # reverse direction
        Z = inverse(odom.Zij)
        Xval = getVal(Xi[2])
    elseif Xid == Xi[2].index
        # forward direction
        Z = odom.Zij
        Xval = getVal(Xi[1])
    else
        error("Bad evalPairwise Pose3Pose3")
    end

    return projectParticles(Xval, Z, odom.Cov)
end

# Still limited to linear sampler, then reprojected onto ball -- TODO upgrade manifold sampler
function evalPotential(odom::Pose3Pose3NH,
      Xi::Array{Graphs.ExVertex,1},
      Xid::Int64;
      N::Int64=200)
   #
  # rz,cz = size(odom.Zij)
  Xval = Array{Float64,2}()
  # implicit equation portion -- bi-directional pairwise function made explicit here
  if Xid == Xi[1].index #odom.
      # reverse direction
      Z = inverse(odom.Zij)
      Xval = getVal(Xi[2])
      oldval = getVal(Xi[1])
  elseif Xid == Xi[2].index
      # forward direction
      Z = odom.Zij
      Xval = getVal(Xi[1])
      oldval = getVal(Xi[2])
  else
      error("Bad evalPairwise Pose3Pose3")
  end

  skipdos = rand(odom.ValidHypot, N)
  dos, donts = skipdos .== 2, skipdos .== 1

  projted = Array{Float64}(6, N)
  projted[:,dos] = projectParticles(Xval[:,dos], Z, odom.Cov)

  projted[:,donts] = oldval[:,donts]

  return projted
end










#
