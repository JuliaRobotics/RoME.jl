# Linear array sonar constraints

# These types should be consolidated in some form -- still exploring for good factorization
"""
$(TYPEDEF)
"""
mutable struct reuseLBRA
  wTb::SE3
  M::Array{Float64,2}
  E::Euler
  inp::Vector{Float64}
  outp::Vector{Float64}
  rbe::Vector{Float64}
  reuseLBRA()=new()
  reuseLBRA(::Int)=new(SE3(0),Matrix{Float64}(LinearAlgebra.I, 4,4),Euler(0),ones(4),ones(4), zeros(3))
  reuseLBRA(a,b,c,d,e,f)=new(a,b,c,d,e,f)
end

"""
$(TYPEDEF)
"""
mutable struct LinearRangeBearingElevation <: FunctorPairwiseMinimize
  range::Normal
  bearing::Normal
  elev::Uniform
  reuse::Vector{reuseLBRA}
  LinearRangeBearingElevation() = new()
  LinearRangeBearingElevation( r::Tuple{Float64,Float64}, b::Tuple{Float64,Float64}; elev=Uniform(-0.25133,0.25133)) = new(Normal(r...),Normal(b...),elev, reuseLBRA[reuseLBRA(0) for i in 1:Threads.nthreads()] )
end
function (p::LinearRangeBearingElevation)(
            res::Vector{Float64},
            userdata::FactorMetadata,
            idx::Int,
            meas::Tuple{Array{Float64,2}},
            pose::Array{Float64,2},
            landm::Array{Float64,2}  )
  #
  residualLRBE!(res, meas[1][:,idx], pose[:,idx], landm[:,idx], p.reuse[Threads.threadid()])
  return res[1]
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
  b = atan(bTl[2],	bTl[1])
  el = -atan(bTl[3], bTl[1])
  return [norm(bTl[1:3]); b; el]
end

function ominus!(reuse::reuseLBRA, X::Vector{Float64}, L::Array{Float64})
  copyto!(reuse.wTb.t, X[1:3])
  reuse.E.R, reuse.E.P, reuse.E.Y = X[4], X[5], X[6]
  convert!(reuse.wTb.R, reuse.E)  # costly
  matrix!(reuse.M, reuse.wTb)
  reuse.inp[1:3] = L[1:3]
  reuse.outp[1:4] = reuse.M\reuse.inp  # bTl  # costly
  reuse.rbe[1] = norm(reuse.outp[1:3])
  reuse.rbe[2] = atan(reuse.outp[2],  reuse.outp[1])
  reuse.rbe[3] = -atan(reuse.outp[3], reuse.outp[1])
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
  resid[1] = norm(z - reuse.rbe)
  # resid[:] = z - ominus!(LinearRangeBearingElevation, X, L)
  nothing
end






#
