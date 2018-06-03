# Linear array sonar constraints

# These types should be consolidated in some form -- still exploring for good factorization
type reuseLBRA
  wTb::SE3
  M::Array{Float64,2}
  E::Euler
  inp::Vector{Float64}
  outp::Vector{Float64}
  rbe::Vector{Float64}
  reuseLBRA()=new()
  reuseLBRA(::Int)=new(SE3(0),eye(4),Euler(0),ones(4),ones(4), zeros(3))
  reuseLBRA(a,b,c,d,e,f)=new(a,b,c,d,e,f)
end
type LinearRangeBearingElevation <: FunctorPairwise
  range::Normal
  bearing::Normal
  elev::Uniform
  reuse::reuseLBRA
  LinearRangeBearingElevation() = new()
  LinearRangeBearingElevation( r::Tuple{Float64,Float64}, b::Tuple{Float64,Float64}; elev=Uniform(-0.25133,0.25133)) = new(Normal(r...),Normal(b...),elev, reuseLBRA(0))
end
function (p::LinearRangeBearingElevation)(
            res::Vector{Float64},
            userdata::Union{Void, FactorMetadata},
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            pose::Array{Float64,2},
            landm::Array{Float64,2}  )
  #
  residualLRBE!(res, meas[1][:,idx], pose[:,idx], landm[:,idx], p.reuse)
  nothing
end
function (p::LinearRangeBearingElevation)(
            res::Vector{Float64},
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            pose::Array{Float64,2},
            landm::Array{Float64,2}  )
  #
  p(res, nothing, idx, meas, pose, landm)
end

function getSample!(y::Array{Float64,2}, las::LinearRangeBearingElevation, idx::Int )
  y[1,idx] = rand(las.range)
  y[2,idx] = rand(las.bearing)
  y[3,idx] = rand(las.elev)
  nothing
end
function getSample( las::LinearRangeBearingElevation, N::Int=1 )
  y = zeros(3,N)
  for i in 1:N
    getSample!(y, las, i)
  end
  return (y,)
end



# returns [Range Bearing Elevation] manifold difference between pose X ominus landmark L
function ominus(::Type{LinearRangeBearingElevation}, X::Vector{Float64}, L::Vector{Float64})
  # rangeBearing3(X, L)
  wTb = SE3(X[1:3], Euler(X[4:6]...))
  bTl = matrix(wTb)\[L[1:3];1.0]
  b = atan2(bTl[2],	bTl[1])
  el = -atan2(bTl[3], bTl[1])
  return [norm(bTl[1:3]); b; el]
end

function ominus!(reuse::reuseLBRA, X::Vector{Float64}, L::Array{Float64})
  copy!(reuse.wTb.t, X[1:3])
  reuse.E.R, reuse.E.P, reuse.E.Y = X[4], X[5], X[6]
  convert!(reuse.wTb.R, reuse.E)  # costly
  matrix!(reuse.M, reuse.wTb)
  reuse.inp[1:3] = L[1:3]
  reuse.outp[1:4] = reuse.M\reuse.inp  # bTl  # costly
  reuse.rbe[1] = norm(reuse.outp[1:3])
  reuse.rbe[2] = atan2(reuse.outp[2],	reuse.outp[1])
  reuse.rbe[3] = -atan2(reuse.outp[3], reuse.outp[1]) #-
  nothing
end

# residual should equal zero when system is in balance
# measurement z is measurement vector with [range; bearing; elevation]
# variables are tuple (pose X [dim6], landmark L [dim3])
# function handle follows required parameter list
function residualLRBE!(resid::Vector{Float64}, z::Vector{Float64}, X::Vector{Float64}, L::Array{Float64}, reuse::reuseLBRA)
  # TODO upgrade so the - sign here is used on a manifold too, ominus(z,  ominus(tt, variables...)  )
  # TODO just switch directly to parameterized function
  ominus!(reuse, X, L)
  resid[1:3] = z - reuse.rbe
  # resid[:] = z - ominus!(LinearRangeBearingElevation, X, L)
  nothing
end


# Convolution of conditional to project landmark location from position X (dim6)
# Y (dim3) are projected points
function project!(meas::LinearRangeBearingElevation, pose::Array{Float64,2}, landmark::Array{Float64,2}, idx::Int, ggtemp::Function)
  z = getSample(meas)
  # warn("didson project! to be upgraded") # TODO
  gg = (res, x) -> residualLRBE!(res, z, pose[1:6,idx][:], x, ggtemp.reuse)
  landmark[1:3,idx] = numericRootGenericRandomizedFnc(gg, 3, 3, landmark[1:3,idx][:])
  # landmark[1:3,idx] = numericRootGenericRandomized(residualLRBE!, 3, getSample(meas), pose[1:6,idx][:], landmark[1:3,idx][:]) # ( bearrange3!,
  # x0[1:3,idx] = numericRoot(bearrange3!, getSample(meas), fixed[1:6,idx], x0[1:3,idx]+0.01*randn(3))
  nothing
end

# zsolving for pose given landmark [P | L]
function backprojectRandomized!(meas::LinearRangeBearingElevation,
        landmark::Array{Float64,2},
        pose::Array{Float64,2},
        idx::Int  )
  #
  z = getSample(meas)
  gg = (res, x) -> residualLRBE!(res, z, x, landmark[1:3,idx][:])
  pose[1:6,idx] = numericRootGenericRandomizedFnc(gg, 3, 6, pose[1:6,idx][:])
  nothing
end


function backprojectRandomized!{T}(meas::LinearRangeBearingElevation,
        landmark::Array{Float64,2},
        pose::Array{Float64,2},
        idx::Int,
        gg::T  )
  #
  gg.landmark[1:3] = landmark[1:3,idx]
  gg.pose[1:6] = pose[1:6,idx]
  gg.z[1:3] = getSample(meas)

  # println("doing this")
  fgr = FastGenericRoot{T}(6, 3, gg)
  #initial guess x0
  copy!(fgr.X, pose[1:6,idx])

  numericRootGenericRandomizedFnc!( fgr )
  pose[1:6,idx] = fgr.Y
	nothing
end

function backprojectRandomized!{T}( fgr::FastGenericRoot{T} )
  #
  # fgr.usrfnc.z[1:3] = getSample(fgr.z)
  numericRootGenericRandomizedFnc!( fgr )
  nothing
end



# convenience function for project!
function project(meas::LinearRangeBearingElevation, X::Vector{Float64}, Y::Vector{Float64})
  xx, yy = X', Y'
  project!(meas, xx', yy', 1)
  return yy[1:3]
end


# (p::WrapParam{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
#     residualLRBE!(res, p.z, x, p.landmark, p.reuse)
# (p::WrapParamArray{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
#     residualLRBE!(res, p.z, x, p.landmark[:,p.idx], p.reuse)




function +(arr::Array{Float64,2}, meas::LinearRangeBearingElevation)
  N = size(arr,2)
  L = zeros(3,N);
  t = Array{Array{Float64,2},1}()
  push!(t,arr)
  push!(t,L)

  fp! = GenericWrapParam{LinearRangeBearingElevation}(meas, t, 2, 1, (zeros(0,1),) , getSample)
  # pre-emptively populate the measurements, kept separate since nlsolve calls fp(x, res) multiple times
  fp!.measurement = fp!.samplerfnc(fp!.usrfnc!, N)
  # fp!(x, res)
  fr = FastRootGenericWrapParam{LinearRangeBearingElevation}(fp!.params[fp!.varidx], 3, fp!)
  for fp!.particleidx in 1:N
    numericRootGenericRandomizedFnc!( fr )
  end
  return L
end
