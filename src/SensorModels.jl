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
mutable struct LinearRangeBearingElevation <: IIF.AbstractManifoldMinimize
  range::Normal
  bearing::Normal
  elev::Uniform
  reuse::Vector{reuseLBRA}
  LinearRangeBearingElevation() = new()
  LinearRangeBearingElevation( r::Tuple{Float64,Float64}, b::Tuple{Float64,Float64}; elev=Uniform(-0.25133,0.25133)) = new(Normal(r...),Normal(b...),elev, reuseLBRA[reuseLBRA(0) for i in 1:Threads.nthreads()] )
end

getManifold(::LinearRangeBearingElevation) = Euclidean(3)#FIXME

function (cfo::CalcFactor{<:LinearRangeBearingElevation})(meas, _pose, _landm)
  #FIXME update to manifolds, quick fix convert for now
  pose = getCoordinates(Pose3, _pose)
  landm = getCoordinates(Point3, _landm)
  return residualLRBE!(meas, pose, landm, cfo.factor.reuse[Threads.threadid()])  
end

function getSample!(y::AbstractVector{<:Real}, las::LinearRangeBearingElevation)
  y[1] = rand(las.range)
  y[2] = rand(las.bearing)
  y[3] = rand(las.elev)
  nothing
end
function getSample( cfo::CalcFactor{<:LinearRangeBearingElevation})
  y = zeros(3)
  getSample!(y, cfo.factor)
  return (y,)
end



# returns [Range Bearing Elevation] manifold difference between pose X ominus landmark L
function ominus(::Type{LinearRangeBearingElevation}, 
                X::AbstractVector{<:Real}, 
                L::AbstractVector{<:Real})
  # rangeBearing3(X, L)
  wTb = SE3(X[1:3], Euler(X[4:6]...))
  bTl = matrix(wTb)\[L[1:3];1.0]
  b = atan(bTl[2],	bTl[1])
  el = -atan(bTl[3], bTl[1])
  return [norm(bTl[1:3]); b; el]
end

function ominus!( reuse::reuseLBRA, 
                  X::AbstractVector{<:Real}, 
                  L::AbstractVector{<:Real} )
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
function residualLRBE!( z::AbstractVector{<:Real}, 
                        X::AbstractVector{<:Real}, 
                        L::AbstractVector{<:Real}, 
                        reuse::reuseLBRA )
  #
  # TODO upgrade so the - sign here is used on a manifold too, ominus(z,  ominus(tt, variables...)  )
  # TODO just switch directly to parameterized function
  
  ominus!(reuse, X, L)
  return  z .- reuse.rbe
end






#
