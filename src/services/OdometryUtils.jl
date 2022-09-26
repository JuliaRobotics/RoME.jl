# utility functions that provide Gaussian odometry accumulation



"""
    $SIGNATURES

Advance an odometry factor as though integrating an ODE -- i.e. ``X_2 = X_1 ⊕ ΔX``. Accepts continuous domain process noise density `Qc` which is internally integrated to discrete process noise Qd.  ``DX`` is assumed to already be incrementally integrated before this function.  See related `accumulateContinuousLocalFrame!` for fully continuous system propagation.

Notes
- This update stays in the same reference frame but updates the local vector as though accumulating measurement values over time.
- Kalman filter would have used for noise propagation: ``Pk1 = F*Pk*F' + Qdk``
- From Chirikjian, Vol.II, 2012, p.35: Jacobian SE(2), Jr = [cθ sθ 0; -sθ cθ 0; 0 0 1] -- i.e. dSE2/dX' = SE2([0;0;-θ])
- `DX = dX/dt*Dt`
- assumed process noise for `{}^b Qc = {}^b [x;y;yaw] = [fwd; sideways; rotation.rate]`

Dev Notes
- TODO many operations here can be done in-place.

Related

accumulateContinuousLocalFrame!, accumulateDiscreteReferenceFrame!, [`accumulateFactorMeans`](@ref)
"""
function accumulateDiscreteLocalFrame!( mpp::MutablePose2Pose2Gaussian,
                                        DX::Vector{Float64},
                                        Qc::Matrix{Float64},
                                        dt::Float64=1.0;
                                        Fk = SE2([0;0;-DX[3]]),
                                        Gk = Matrix{Float64}(LinearAlgebra.I, 3,3),
                                        Phik = SE2(DX) )
  #
  kXk1 = SE2(mpp.Z.μ)*Phik
  phi, gamma, Qd = cont2disc(Fk, Gk, Qc, dt, Phik)
  Covk1 = Phik*(mpp.Z.Σ.mat)*(Phik') + Qd
  check = norm(Covk1 - Covk1')
  1e-4 < check ? @warn("Covk1 is getting dangerously non-Hermitian, still forcing symmetric covariance matrix.") : nothing
  @assert check < 1.0
  Covk1 .+= Covk1'
  Covk1 ./= 2
  mpp.Z = MvNormal(se2vee(kXk1), Covk1)
  nothing
end

accumulateDiscreteLocalFrame!(dfg::AbstractDFG,
                              fctlbl::Symbol,
                              DX::Vector{Float64},
                              Qc::Matrix{Float64},
                              dt::Float64=1.0;
                              Fk = SE2([0;0;-DX[3]]),
                              Gk = Matrix{Float64}(LinearAlgebra.I, 3,3),
                              Phik = SE2(DX) ) = accumulateDiscreteLocalFrame!(getFactorFunction(dfg, fctlbl),DX, Qc, dt; Fk=Fk, Gk=Gk, Phik=Phik)

"""
    $SIGNATURES

Helper function to duplicate values from a special factor variable into standard factor and variable.  Returns the name of the new factor.

Notes:
- Developed for accumulating odometry in a `MutablePosePose` and then cloning out a standard PosePose and new variable.
- Does not change the original MutablePosePose source factor or variable in any way.
- Assumes timestampe from mpp object.

Related

[`addVariable!`](@ref), [`addFactor!`](@ref)
"""
function duplicateToStandardFactorVariable( ::Type{Pose2Pose2},
                                            mpp::MutablePose2Pose2Gaussian,
                                            dfg::AbstractDFG,
                                            prevsym::Symbol,
                                            newsym::Symbol;
                                            solvable::Int=1,
                                            graphinit::Bool=true,
                                            cov::Union{Nothing, Matrix{Float64}}=nothing  )::Symbol
  #

  # extract factor values and create PosePose object
  posepose = Pose2Pose2(MvNormal(mpp.Z.μ, cov===nothing ? mpp.Z.Σ.mat : cov))

  # modify the factor graph
  addVariable!(dfg, newsym, Pose2, solvable=solvable, timestamp=mpp.timestamp)
  fct = addFactor!(dfg, [prevsym; newsym], posepose, solvable=solvable, graphinit=graphinit, timestamp=mpp.timestamp)
  # new factor name
  # return ls(dfg, newsym)[1]
  return DFG.getLabel(fct)
end

"""
    $SIGNATURES

Reset the transform value stored in a `::MutablePose2Pose2Gaussian` to zero.
"""
function resetFactor!(mpp::MutablePose2Pose2Gaussian)::Nothing
  mpp.Z = MvNormal(zeros(3), 1e-6*Matrix{Float64}(LinearAlgebra.I, 3,3) )
  nothing
end


"""
    $SIGNATURES

Extract deltas from existing dead reckoning data so that odometry calculations can be repeated later.

Notes
- Useful for reverse engineering data or simulation tools.
"""
function extractDeltaOdo(XX, YY, TH)
  dt = 1.0
  DX = zeros(3,length(XX))
  nXYT__ = zeros(3,size(DX,2))
  nXYT__[:,1] = [XX[1];YY[1];TH[1]]
  for i in 2:length(XX)
    wTbk = SE2([XX[i-1];YY[i-1];TH[i-1]])
    wTbk1 = SE2([XX[i];YY[i];TH[i]])
    bkTbk1 = wTbk\wTbk1
    DX[:,i] = se2vee(bkTbk1)

    # test
    nXYT__[:,i] .= se2vee(SE2(nXYT__[:,i-1])*SE2(DX[:,i]))
    # nXYT__[:,i] .= se2vee(SE2(nXYT__[:,i-1])*bkTbk1)
  end

  return DX
end


# DX = [transx, transy, theta]
function addPose2Pose2!(retval::Array{<:Real,1}, x::Array{<:Real,1}, dx::Array{<:Real,1})
  X = SE2(x)
  DX = SE2(dx)
  se2vee!(retval, X*DX)
  nothing
end
function addPose2Pose2(x::Array{<:Real,1}, dx::Array{<:Real,1})
  retval = zeros(3)
  addPose2Pose2!(retval, x, dx)
  return retval
end


## Previous methods

function odomKDE(p1,dx,cov)
  @warn "odomKDE is beig deprecated in its current form, consider using approxConv or predictVariableByFactor instead."
  X = getPoints(p1)
  sig = diag(cov)
  RES = zeros(size(X))
  # increases the number of particles based on the number of modes in the measurement Z
  for i in 1:size(X,2)
      ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
      RES[:,i] = addPose2Pose2(X[:,i], dx + ent)
  end
  return manikde!(RES, Pose2)
end

"""
    $SIGNATURES

Calculate the relative chords between consecutive poses in the factor graph.
Data structure is Dict{Symbol,Dict{Symbol,Tuple{Matrix,Matrix}}}.
The two Matrix values are 3x100, with the first as shown in the attached screen capture.
These values should be the relative transform from dict[:x0][:x1], or dict[:x0][:x2], or dict[:x0][:x2] etc for all poses up to some reasonable chord length.
There are also two matrix values: the first is the relative transform based on measurements only, the second matrix is the same relative transform but according to the SLAM solution of any and all data being used.
"""
function assembleChordsDict(dfg::AbstractDFG,
                            vsyms = ls(dfg, r"x\d") |> sortDFG;
                            MAXADI = 10,
                            lastPoseNum = getVariableLabelNumber(vsyms[end]),
                            chords = Dict{Symbol,Dict{Symbol,Tuple}}()  )
  #
  # fsyms = [:x0x1f1; :x1x2f1]


  @sync for from in vsyms[1:end-1]
    SRT = getVariableLabelNumber(from)
    chords[from] = Dict{Symbol,Tuple}()
    maxadi = lastPoseNum - getVariableLabelNumber(from)
    maxadi = MAXADI < maxadi ? MAXADI : maxadi
    for adi in 1:maxadi
      to = Symbol("x",getVariableLabelNumber(from)+adi)
      # FIXME replace with approxConvBelief instead
      tt = Threads.@spawn accumulateFactorChain(dfg, $from, $to)
      @async begin
        chords[$from][$to] = fetch(tt)
      end
    end
  end

  chords
end




"""
    $(SIGNATURES)

Create a new variable node and insert odometry constraint factor between
which will automatically increment latest pose symbol x<k+1> for new node new node and
constraint factor are returned as a tuple.
"""
function addOdoFG!(
  fg::AbstractDFG,
  n::Symbol,
  DX::Array{Float64,1},
  cov::Array{Float64,2};
  N::Int=0,
  solvable::Int=1,
  labels::Vector{<:AbstractString}=String[]  
)
  #
  prev, X, nextn = getLastPose2D(fg)
  r,c = size(X)
  if N==0
    N = c
  end
  sig = diag(cov)
  XnextInit = zeros(r,c)
  # increases the number of particles based on the number of modes in the measurement Z
  for i in 1:c
    ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
    XnextInit[:,i] = addPose2Pose2(X[:,i], DX + ent)
  end

  v = addVariable!(fg, n, Pose2, N=N, solvable=solvable, tags=[labels;"POSE"])
  # v = addVariable!(fg, n, XnextInit, cov, N=N, solvable=solvable, tags=labels)
  pp = Pose2Pose2(MvNormal(DX, cov)) #[prev;v],
  f = addFactor!(fg, [prev;v], pp, solvable=solvable, graphinit=true )
  infor = inv(cov^2)
  # addOdoRemote(prev.index,v.index,DX,infor) # this is for remote factor graph ref parametric solution -- skipped internally by global flag variable
  return v, f
end

function addOdoFG!(
  fgl::AbstractDFG,
  Z::Pose3Pose3;
  N::Int=0,
  solvable::Int=1,
  labels::Vector{<:AbstractString}=String[]  
)
  #
  vprev, X, nextn = getLastPoses(fgl)[1]
  vnext = addVariable!(fgl, nextn, Pose3, solvable=solvable, tags=labels)
  fact = addFactor!(fgl, [vprev;vnext], Z, graphinit=true)

  return vnext, fact
  # error("addOdoFG!( , ::Pose3Pose3, ) not currently usable, there were breaking changes. Work in Progress")
  # addOdoFG(fg, n, DX, cov, N=N, solvable=solvable, tags=labels)
end

"""
    $(SIGNATURES)

Create a new variable node and insert odometry constraint factor between
which will automatically increment latest pose symbol x<k+1> for new node new node and
constraint factor are returned as a tuple.

"""
function addOdoFG!(
  fgl::AbstractDFG,
  odo::Pose2Pose2;
  N::Int=0,
  solvable::Int=1,
  labels::Vector{<:AbstractString}=String[] 
)
  #
  vprev, X, nextn = getLastPose(fgl)
  if N==0
    N = size(X,2)
  end
  # vnext = addVariable!(fgl, nextn, X⊕odo, ones(1,1), N=N, solvable=solvable, tags=labels)
  vnext = addVariable!(fgl, nextn, Pose2, N=N, solvable=solvable, tags=labels)
  fact = addFactor!(fgl, [vprev;vnext], odo, graphinit=true)

  return vnext, fact
end


#
