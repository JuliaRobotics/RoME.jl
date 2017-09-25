
# And 2xN matrix of N landmark xy positions as variable nodes in factor graph
function addSimMapFG!(fg::FactorGraph, lms::Array{Float64,2})
    for i in 1:size(lms,2)
      newLandm!(fg,string('l',i), vectoarr2(lms[:,i]), 0.001*eye(2))
    end
    nothing
end

# assume all the landmarks are already loaded into the Ground Truth FG
function simOdo!(fgGT::FactorGraph, fg::FactorGraph, DX::Array{Float64,1};
    noiserate=2.0*[3e-2;3e-2;1.5e-3], driftrate=[0.0;0.0;0.0], detLM=Union{})
    prev, X, nextn = getLastPose2D(fg)
    addOdoFG!(fgGT, nextn, DX, 0.001*eye(3))

    r = norm(DX[1:2])
    xn = noiserate[1]*r
    yn = noiserate[2]*r
    thn = noiserate[3]*r
    cov = diagm([xn;yn;thn])
    DXn = DX + [xn*randn();yn*randn();thn*randn()] + r*driftrate
    addOdoFG!(fg, nextn, DXn, cov)

    return nextn
end


function truePredBR(fgGT::FactorGraph, fg::FactorGraph, ps::String, lm::String)
    trubr = predictBodyBR(fgGT, ps, lm)
    truA = [trubr[1]; trubr[2]]
    # Prediction of BR measurement
    prebr = predictBodyBR(fg, ps, lm)
    preA = [prebr[1];prebr[2]]
    return truA, preA
end

function showTruePredBR(fgGT::FactorGraph, fg::FactorGraph, ps::String, lm::String, cov::Array{Float64,2})
    truA, preA = truePredBR(fgGT, fg, ps, lm)
    measA = truA + [cov[1,1]*randn();cov[2,2]*randn()]
    mala = malahanobisBR(measA, preA, cov)
    println(ps,lm, ": true BR=$(round(truA,3)), pred BR=$(round(preA,3)), mala=$(round(mala,3))")
end

function crossMalaBR(fgGT::FactorGraph, fg::FactorGraph,
                      ps::String, lmT::String, lmE, cov::Array{Float64,2})
    truT, preT = truePredBR(fgGT, fg, ps, lmT)
    measT = truT + [cov[1,1]*randn();cov[2,2]*randn()]

    truE, preE = truePredBR(fgGT, fg, ps, lmE)
    malaX = malahanobisBR(measT, preE, cov)
    return malaX
end
