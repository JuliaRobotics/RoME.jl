

function predictBodyBR(fg::AbstractDFG, a::Symbol, b::Symbol)
  res = zeros(2)
  A = getKDEMean(getBelief(getVariable(fg,a)))
  B = getKDEMean(getBelief(getVariable(fg,b)))
  Ax = A[1] # Statistics.mean(vec(A[1,:]))
  Ay = A[2] # Statistics.mean(vec(A[2,:]))
  Ath = getKDEMax(getBelief(getVariable(fg,a)))[3]
  Bx = B[1]
  By = B[2]
  wL = SE2([Bx;By;0.0])
  wBb = SE2([Ax;Ay;Ath])
  bL = se2vee((wBb \ Matrix{Float64}(LinearAlgebra.I, 3,3)) * wL)
  dx = bL[1] - 0.0
  dy = bL[2] - 0.0
  b = (atan(dy,dx))
  r = sqrt(dx^2 + dy^2)
  return b, r
end


function malahanobisBR(measA, preA, cov::Array{Float64,2})
  # measure landmark with noise
  res = measA - preA
  mala2 = Union{}
  #Malahanobis distance
  if false
    lambda = cov \ Matrix{Float64}(LinearAlgebra.I, 2,2)
    mala2 = res' * lambda * res
  else
    mala2 = res' * (cov \ res)
  end

  mala = sqrt(mala2)
  return mala
end



"""
    $SIGNATURES

Method to compare current and predicted estimate on a variable, developed for testing a new factor before adding to the factor graph.

Notes
- `fct` does not have to be in the factor graph -- likely used to test beforehand.
- function is useful for detecting if `multihypo` should be used.
- `approxConv` will project the full belief estimate through some factor but must already be in factor graph.

Example

```julia
# fg already exists containing :x7 and :l3
pp = Pose2Point2BearingRange(Normal(0,0.1),Normal(10,1.0))
# possible new measurement from :x7 to :l3
curr, pred = predictVariableByFactor(fg, :l3, pp, [:x7; :l3])
# example of naive user defined test on fit score
fitscore = minkld(curr, pred)
# `multihypo` can be used as option between existing or new variables
```

Related

approxConv
"""
function predictVariableByFactor( dfg::AbstractDFG,
                                  targetsym::Symbol,
                                  fct::Pose2Point2BearingRange,
                                  prevars::Vector{Symbol}  )
  #
  @assert targetsym in prevars
  curr = getBelief(dfg, targetsym)
  tfg = initfg()
  for var in prevars
    varnode = getVariable(dfg, var)
    addVariable!(tfg, var, getSofttype(varnode))
    if var != targetsym
      @assert isInitialized(varnode)
      initVariable!(tfg,var,getBelief(varnode))
    end
  end
  addFactor!(tfg, prevars, fct, graphinit=false)
  fctsym = ls(tfg, targetsym)

  pts, infd = predictbelief(tfg, targetsym, fctsym)
  pred = manikde!(getManifold(getVariable(dfg, targetsym)), pts)
  # return current and predicted beliefs
  return curr, pred
end
