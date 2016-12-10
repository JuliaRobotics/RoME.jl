# Pose2 functions for Robot Motion Estimate


type PriorPose2 <: IncrementalInference.Singleton
    Zi::Array{Float64,2}
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    PriorPose2() = new()
    PriorPose2(x...) = new(x[1], x[2], x[3])
end

type Pose2Pose2 <: IncrementalInference.Pairwise
    Zij::Array{Float64,2} # 2translations, 1rotation
    Cov::Array{Float64,2}
    W::Array{Float64,1}
    Pose2Pose2() = new()
    Pose2Pose2(x...) = new(x[1], x[2], x[3])
end



function evalPotential(obs::PriorPose2, Xi::Array{Graphs.ExVertex,1}; N::Int64=200)
    cov = diag(obs.Cov)
    ret = zeros(3,N)
    for j in 1:N
      for i in 1:size(obs.Zi,1)
        ret[i,j] += obs.Zi[i,1] + (cov[i]*randn())
      end
    end
    return ret
end


# TODO Switch to using SE(2) oplus
# DX = [transx, transy, theta]
function addPose2Pose2(x::Array{Float64,1}, dx::Array{Float64,1})
    X = SE2(x)
    DX = SE2(dx)
    return se2vee(X*DX)
end

function addPose2Pose2!(retval::Array{Float64,1}, x::Array{Float64,1}, dx::Array{Float64,1})
  X = SE2(x)
  DX = SE2(dx)
  se2vee!(retval, X*DX)
  nothing
end

function evalPotential(odom::Pose2Pose2, Xi::Array{Graphs.ExVertex,1}, Xid::Int64)
    rz,cz = size(odom.Zij)
    Xval = Array{Float64,2}()
    XvalNull = Array{Float64,2}()
    # implicit equation portion -- bi-directional pairwise function
    if Xid == Xi[1].index #odom.
        #Z = (odom.Zij\eye(rz)) # this will be used for group operations
        Z = se2vee(SE2(vec(odom.Zij)) \ eye(3))
        Xval = getVal(Xi[2])
        XvalNull = getVal(Xi[1])
    elseif Xid == Xi[2].index
        Z = odom.Zij
        Xval = getVal(Xi[1])
        XvalNull = getVal(Xi[2])
    else
        error("Bad evalPairwise Pose2Pose2")
    end

    r,c = size(Xval)
    RES = zeros(r,c*cz)

    # TODO -- this should be the covariate error from Distributions, only using diagonals here (approxmition for speed in first implementation)
    # dd = size(Z,1) == r
    ENT = randn(r,c)
    HYP = rand(Categorical(odom.W),c) # TODO consolidate
    HYP -= length(odom.W)>1 ? 1 : 0
    for d in 1:r
       @fastmath @inbounds ENT[d,:] = ENT[d,:].*odom.Cov[d,d]
    end
    # Repeat ENT values for new modes from meas
    for j in 1:cz
      for i in 1:c
        if HYP[i]==1 # TODO consolidate hypotheses on Categorical
          z = Z[1:r,j].+ENT[1:r,i]
          RES[1:r,i*j] = addPose2Pose2(Xval[1:r,i], z )
        else
          RES[1:r,i*j] = XvalNull[1:r,i]
        end
      end
    end

    return RES
end



# Temporaru functions for database support -- will be reduced to macro in future
# ------------------------------------


type PackedPriorPose2  <: IncrementalInference.PackedInferenceType
    vecZij::Array{Float64,1} # 0rotations, 1translation in each column
    dimz::Int64
    vecCov::Array{Float64,1}
    dimc::Int64
    W::Array{Float64,1}
    PackedPriorPose2() = new()
    PackedPriorPose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  Zi = reshapeVec2Mat(d.vecZij,d.dimz)
  Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  return PriorPose2(Zi, Cov, d.W)
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  v1 = d.Zi[:];
  v2 = d.Cov[:];
  return PackedPriorPose2(v1,size(d.Zi,1),
                          v2,size(d.Cov,1),
                          d.W)
end
# function convert(::Type{PackedFunctionNodeData{PackedPriorPose2}}, d::FunctionNodeData{PriorPose2})
#   return PackedFunctionNodeData{PackedPriorPose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           string(d.frommodule), convert(PackedPriorPose2, d.fnc))
# end
# function convert(::Type{FunctionNodeData{PriorPose2}}, d::PackedFunctionNodeData{PackedPriorPose2})
#   return FunctionNodeData{PriorPose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           Symbol(d.frommodule), convert(PriorPose2, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{PriorPose2})
#   return convert(PackedFunctionNodeData{PackedPriorPose2}, d)
# end
# function FNDdecode(d::PackedFunctionNodeData{PackedPriorPose2})
#   return convert(FunctionNodeData{PriorPose2}, d)
# end


# --------------------------------------------


type PackedPose2Pose2  <: IncrementalInference.PackedInferenceType
  vecZij::Array{Float64,1} # 2translations, 1rotation
  dimz::Int64
  vecCov::Array{Float64,1}
  dimc::Int64
  W::Array{Float64,1}
  PackedPose2Pose2() = new()
  PackedPose2Pose2(x...) = new(x[1], x[2], x[3], x[4], x[5])
end
function convert(::Type{Pose2Pose2}, d::PackedPose2Pose2)
  return Pose2Pose2(reshapeVec2Mat(d.vecZij,d.dimz),
                    reshapeVec2Mat(d.vecCov, d.dimc), d.W)
end
function convert(::Type{PackedPose2Pose2}, d::Pose2Pose2)
  v1 = d.Zij[:];
  v2 = d.Cov[:];
  return PackedPose2Pose2(v1,size(d.Zij,1),
                          v2,size(d.Cov,1),
                          d.W)
end
# function convert(::Type{PackedFunctionNodeData{PackedPose2Pose2}}, d::FunctionNodeData{Pose2Pose2})
#   return PackedFunctionNodeData{PackedPose2Pose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           string(d.frommodule), convert(PackedPose2Pose2, d.fnc))
# end
# function convert(::Type{FunctionNodeData{Pose2Pose2}}, d::PackedFunctionNodeData{PackedPose2Pose2})
#   return FunctionNodeData{Pose2Pose2}(d.fncargvID, d.eliminated, d.potentialused, d.edgeIDs,
#           Symbol(d.frommodule), convert(Pose2Pose2, d.fnc))
# end
# function FNDencode(d::FunctionNodeData{Pose2Pose2})
#   return convert(PackedFunctionNodeData{PackedPose2Pose2}, d)
# end
# function FNDdecode(d::PackedFunctionNodeData{PackedPose2Pose2})
#   return convert(FunctionNodeData{Pose2Pose2}, d)
# end


# --------------------------------------------
