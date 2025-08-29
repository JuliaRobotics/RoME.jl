
struct P2P2BearingReuse
    measvec::Vector{Float64}
    predvec::Vector{Float64}
    resid::Vector{Float64}
    P2P2BearingReuse() = new(zeros(2), zeros(2), zeros(2))
end

"""
    $TYPEDEF

Single dimension bearing constraint from Pose2 to Point2 variable.
"""
DFG.@defObservationType Pose2Point2Bearing RelativeObservation SpecialOrthogonalGroup(2)

function preambleCache(
    ::AbstractDFG,
    ::AbstractVector{<:VariableCompute},
    ::Pose2Point2Bearing,
)
    return P2P2BearingReuse()
end

function (cfo::CalcFactor{<:Pose2Point2Bearing})(X, p, l)
    # wl = l
    # wTp = p
    # pl = pTw*wl

    pl = transpose(p.x[2]) * (l - p.x[1])
    # δθ = mθ - plθ  # X[2,1] because we store [measurement tangent as 2x2 matrix (Lie algebra, so(2))](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/group.html#Manifolds.exp_lie-Tuple{Manifolds.GeneralUnitaryMultiplicationGroup{ManifoldsBase.TypeParameter{Tuple{2}},%20%E2%84%9D},%20Any})
    δθ = Manifolds.sym_rem(X[2, 1] - atan(pl[2], pl[1]))
    return [δθ]
end
