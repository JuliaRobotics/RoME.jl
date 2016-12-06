# Pose3Pose3 evaluation functions





type PriorPose3 <: IncrementalInference.Singleton
    Zi::SE3
    Cov::Array{Float64,2}
    PriorPose3() = new()
    PriorPose3(st::FloatInt, sr::Float64) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]])
    PriorPose3(s::SE3, c::Array{Float64,2}) = new(s,c)
end
# type PackedPriorPose3
#     vecZij::Array{Float64,1} # 0rotations, 1translation in each column
#     dimz::Int64
#     vecCov::Array{Float64,1}
#     dimc::Int64
#     W::Array{Float64,1}
#     PackedPriorPose3() = new()
#     PackedPriorPose3(x...) = new(x[1], x[2], x[3], x[4], x[5])
# end
# function convert(::Type{PriorPose3}, d::PackedPriorPose3)
#   Zij = reshapeVec2Mat(d.vecZij,d.dimz)
#   Cov = reshapeVec2Mat(d.vecCov, d.dimc)
#   return PriorPose3(Zij, Cov, d.W)
# end
# function convert(::Type{PackedPriorPose3}, d::PriorPose3)
#   v1 = d.Zij[:];
#   v2 = d.Cov[:];
#   return PackedPriorPose3(v1,size(d.Zij,1),
#                           v2,size(d.Cov,1),
#                           d.W)
# end
# function convert(::Type{FunctionNodeData{PackedPriorPose3}}, d::FunctionNodeData{PriorPose3})
#   return FunctionNodeData{PackedPriorPose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PackedPriorPose3, d.fnc))
# end
# function convert(::Type{FunctionNodeData{PriorPose3}}, d::FunctionNodeData{PackedPriorPose3})
#   return FunctionNodeData{PriorPose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PriorPose3, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{PriorPose3})
#   return convert(FunctionNodeData{PackedPriorPose3}, d)
# end
# function FNDdecode(d::FunctionNodeData{PackedPriorPose3})
#   return convert(FunctionNodeData{PriorPose3}, d)
# end



# ------------------------------------


type Pose3Pose3 <: IncrementalInference.Pairwise
    Zij::SE3 # 3translations, 3exponential param rotation
    Cov::Array{Float64,2}
    Pose3Pose3() = new()
    Pose3Pose3(st::FloatInt, sr::Float64) = new(SE3(0), [[st*eye(3);zeros(3,3)];[zeros(3);sr*eye(3)]])
    Pose3Pose3(s::SE3, c::Array{Float64,2}) = new(s,c)
end
# type PackedPose3Pose3
#   vecZij::Array{Float64,1} # 2translations, 1rotation
#   dimz::Int64
#   vecCov::Array{Float64,1}
#   dimc::Int64
#   W::Array{Float64,1}
#   PackedPose2Pose2() = new()
#   PackedPose2Pose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
# end
# function convert(::Type{Pose3Pose3}, d::PackedPose3Pose3)
#   return Pose3Pose3(reshapeVec2Mat(d.vecZij,d.dimz),
#                     reshapeVec2Mat(d.vecCov, d.dimc), W)
# end
# function convert(::Type{PackedPose3Pose3}, d::Pose3Pose3)
#   v1 = d.Zij[:];
#   v2 = d.Cov[:];
#   return PackedPose3Pose3(v1,size(d.Zij,1),
#                           v2,size(d.Cov,1),
#                           d.W)
# end
# function convert(::Type{FunctionNodeData{PackedPose3Pose3}}, d::FunctionNodeData{Pose3Pose3})
#   return FunctionNodeData{PackedPose3Pose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(PackedPose3Pose3, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Pose3Pose3}}, d::FunctionNodeData{PackedPose3Pose3})
#   return FunctionNodeData{Pose3Pose3}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           convert(Pose3Pose3, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{Pose3Pose3})
#   return convert(FunctionNodeData{PackedPose3Pose3}, d)
# end
# function FNDdecode(d::FunctionNodeData{PackedPose3Pose3})
#   return convert(FunctionNodeData{Pose3Pose3}, d)
# end



# -----------------------






# 3 dimensionsal evaluation functions.
# Transforms code is scattered, and will be moved to better location with stable 3D Examples
# also worth it to consolidate if other affine transforms package is available

# using Euler angles for linear sampling in product of potentials
# Gaussian model for prior
function evalPotential(obs::PriorPose3, Xi::Array{Graphs.ExVertex,1}; N::Int64=200)
  mu = veeEuler(obs.Zi)
  return rand( MvNormal(mu, obs.Cov), N)
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
      x.R, x.t = TransformUtils.convert(SO3,Euler((Xval[4:6,i][:])...)), Xval[1:3,i][:]
      ent.R, ent.t = TransformUtils.convert(SO3, so3(ENT[4:6,i][:])), ENT[1:3,i][:]
      newval = x*Z*ent
      RES[1:r,i*j] = veeEuler(newval)
    end
  # end
  #
  return RES
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
