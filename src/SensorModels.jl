# Linear array sonar constraints

# typealias FloatInt Union{Float64, Int64}

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
      idx::Int,
      meas::Tuple{Array{Float64,2}},
      pose::Array{Float64,2},
      landm::Array{Float64,2}  )
  #
  residualLRBE!(res, meas[1][:,idx], pose[:,idx], landm[:,idx], p.reuse)
  nothing
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

warn("Obsolete definitions of WrapParam and WrapParamArray")
type WrapParam{T} <: Function
  # params::Tuple
  landmark::Vector{Float64}
  pose::Vector{Float64}
  z::Vector{Float64}
  reuse::T
end
# towards Tuple of pointers for speed and flexibility
type WrapParamArray{T} <: Function
  # params::Tuple
  landmark::Array{Float64}
  pose::Array{Float64}
  z::Vector{Float64}
  idx::Int
  reuse::T
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
  # return [norm(reuse.outp[1:3]); b; el]
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
# (p::WrapParam{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) = residualLRBE!(res, p.z, x, p.landmark, p.reuse)
# function (p::LinearRangeBearingElevation)(
#       res::Vector{Float64},
#       )
#   #
#   residualLRBE!(resid, )
#   res[1:3] = p.reuse.resid[1:3]
#   nothing
# end


# Convolution of conditional to project landmark location from position X (dim6)
# Y (dim3) are projected points
function project!(meas::LinearRangeBearingElevation, pose::Array{Float64,2}, landmark::Array{Float64,2}, idx::Int, ggtemp::Function)
  z = getSample(meas)
  # warn("didson project! to be upgraded") # TODO
  gg = (x, res) -> residualLRBE!(res, z, pose[1:6,idx][:], x, ggtemp.reuse)
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
  # error("outdated")
  z = getSample(meas)
  gg = (x, res) -> residualLRBE!(res, z, x, landmark[1:3,idx][:])
  pose[1:6,idx] = numericRootGenericRandomizedFnc(gg, 3, 6, pose[1:6,idx][:])
	# pose[1:6,idx] = numericRootGenericRandomized(residualLRBE!, 3, getSample(meas), landmark[1:3,idx][:], pose[1:6,idx][:]) # ( bearrange3!,
	nothing
end

function backprojectRandomizedOld!(meas::LinearRangeBearingElevation,
        landmark::Array{Float64,2},
        pose::Array{Float64,2},
        idx::Int,
        gg::Function  )
  #
  gg.landmark[1:3] = landmark[1:3,idx][:]
  gg.pose[1:6] = pose[1:6,idx][:]
  gg.z[1:3] = getSample(meas)

  # gg = (x, res) -> residualLRBE!(res, z, x, landmark[1:3,idx][:])
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


(p::WrapParam{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
    residualLRBE!(res, p.z, x, p.landmark, p.reuse)
(p::WrapParamArray{reuseLBRA})(x::Vector{Float64}, res::Vector{Float64}) =
    residualLRBE!(res, p.z, x, p.landmark[:,p.idx], p.reuse)




function evalPotential(meas::LinearRangeBearingElevation, Xi::Array{Graphs.ExVertex,1}, Xid::Int64; N::Int=100)
  # Function soon to be replaced with improved functor version
  fromX, ret, ff = zeros(0,0), zeros(0,0), +

  fp! = WrapParam{reuseLBRA}(zeros(3), zeros(6), zeros(3), reuseLBRA(0))

  if Xi[1].index == Xid
    fromX = getVal( Xi[2] )
    ret = deepcopy(getVal( Xi[1] ))
    ff = backprojectRandomized!
  elseif Xi[2].index == Xid
    fromX = getVal( Xi[1] )
    ret = deepcopy(getVal( Xi[2] ))
    ff = project!
  end
  r,c = size(fromX)

	for i in 1:200
		ff(meas, fromX, ret, i, fp!)
	end

  return ret
end


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

## old code

# function backprojectRandomized!(meas::LinearRangeBearingElevation, landmark::Array{Float64,2}, pose::Array{Float64,2}, idx::Int)
# .....
# Zbr = sample(meas)
#
#   p = collect(1:6);
#   shuffle!(p);
#   p1 = p.==1; p2 = p.==2; p3 = p.==3
#   #@show x0, par
#   r = nlsolve(    (x, res) -> bearrangDidson!(res, Zbr,
#                   shuffleXAltD(x, x0[:,idx][:], 3, p), fixed[1:3,idx][:] ),
#                   [x0[p1,idx];x0[p2,idx];x0[p3,idx]]+0.01*randn(3)   )
#
#   X[1:6,idx] = shuffleXAltD(r.zero, x0[:,idx][:], 3, p );

# # randomized backproject using residual function
# function bearrange3!(residual::Vector{Float64}, Zrb::Vector{Float64}, X::Vector{Float64}, L::Vector{Float64})
#   wTb = SE3(X[1:3], Euler(X[4:6]...))
#   bTl = matrix(wTb)\[L[1:3];1.0]
#   b = atan2(bTl[2],	bTl[1])
# 	el = -atan2(bTl[3], bTl[1])
#   residual[1] = Zrb[1]-norm(bTl[1:3])
#   residual[2] = Zrb[2]-b
# 	residual[3] = Zrb[3]-el
#   nothing
# end

# returns vector with [range, bearing, elevation]
# takes pose X [dim6] and landmark L [dim3]
# function rangeBearing3( X::Vector{Float64}, L::Vector{Float64})
#   wTb = SE3(X[1:3], Euler(X[4:6]...))
#   bTl = matrix(wTb)\[L[1:3];1.0]
#   b = atan2(bTl[2],	bTl[1])
# 	el = -atan2(bTl[3], bTl[1])
#   return [norm(bTl[1:3]); b; el]
# end


# function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
# 		n = length(X)
#     Y = deepcopy(Alt)
# 		for i in 1:d
# 			Y[p[i]] = X[i]
# 		end
#     return Y
# end
#





# residual function must have the form
# bearrangDidson!(residual::Array{Float64,1}, Zbr::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
# function solveLandmDidsonPt(measurement::Array{Float64,1}, param::Array{Float64,1}, x0::Array{Float64,1})
# 	numericRoot(bearrangDidson!, measurement, param, x0)
#   # return (nlsolve(   (X, residual) -> bearrangDidson!(residual, measurement, param, X), x0 )).zero
# end
