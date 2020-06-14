
# And 2xN matrix of N landmark xy positions as variable nodes in factor graph
function addSimMapFG!(fg::AbstractDFG, lms::Array{Float64,2})
    for i in 1:size(lms,2)
      newLandm!(fg,string('l',i), vectoarr2(lms[:,i]), 0.001*Matrix{Float64}(LinearAlgebra.I, 2,2))
    end
    nothing
end

# assume all the landmarks are already loaded into the Ground Truth FG
function simOdo!(fgGT::AbstractDFG, fg::AbstractDFG, DX::Array{Float64,1};
    noiserate=2.0*[3e-2;3e-2;1.5e-3], driftrate=[0.0;0.0;0.0], detLM=Union{})
    prev, X, nextn = getLastPose2D(fg)
    addOdoFG!(fgGT, nextn, DX, 0.001*Matrix{Float64}(LinearAlgebra.I, 3,3))

    r = norm(DX[1:2])
    xn = noiserate[1]*r
    yn = noiserate[2]*r
    thn = noiserate[3]*r
    cov = Matrix(Diagonal([xn;yn;thn]))
    DXn = DX + [xn*randn();yn*randn();thn*randn()] + r*driftrate
    addOdoFG!(fg, nextn, DXn, cov)

    return nextn
end


function truePredBR(fgGT::AbstractDFG, fg::AbstractDFG, ps::String, lm::String)
    trubr = predictBodyBR(fgGT, ps, lm)
    truA = [trubr[1]; trubr[2]]
    # Prediction of BR measurement
    prebr = predictBodyBR(fg, ps, lm)
    preA = [prebr[1];prebr[2]]
    return truA, preA
end

"""
    $SIGNATURES

Calculate a bearing and range parameter between Pose2 and Point2 `::Vector`s.

Example
```julia
b,r = calcPosePointBearingRange([0;0;0], [10;10])
```
"""
function calcPosePointBearingRange(pose::Vector{<:Real},
                                   poin::Vector{<:Real}  )
  #
  @assert length(pose) == 3
  @assert length(poin) == 2

  #
  DD = poin-pose[1:2]
  ran = norm(DD)
  phi = atan(DD[2], DD[1])
  the = TU.wrapRad(phi - pose[3])

  #
  return the, ran
end


function showTruePredBR(fgGT::AbstractDFG, fg::AbstractDFG, ps::String, lm::String, cov::Array{Float64,2})
    truA, preA = truePredBR(fgGT, fg, ps, lm)
    measA = truA + [cov[1,1]*randn();cov[2,2]*randn()]
    mala = malahanobisBR(measA, preA, cov)
    println(ps,lm, ": true BR=$(round(truA,digits=3)), pred BR=$(round(preA,digits=3)), mala=$(round(mala,digits=3))")
end

function crossMalaBR(fgGT::AbstractDFG, fg::AbstractDFG,
                      ps::String, lmT::String, lmE, cov::Array{Float64,2})
    truT, preT = truePredBR(fgGT, fg, ps, lmT)
    measT = truT + [cov[1,1]*randn();cov[2,2]*randn()]

    truE, preE = truePredBR(fgGT, fg, ps, lmE)
    malaX = malahanobisBR(measT, preE, cov)
    return malaX
end
