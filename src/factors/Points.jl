
(cf::CalcFactor{<:PriorPoint2})(X, p) = X[1:2] .- p[1:2]

(pp2r::CalcFactor{<:Point2Point2})(X, p, q) = X[1:2] .- (q[1:2] .- p[1:2])

function (cf::CalcFactor{<:PriorPoint3})(X, p::ArrayPartition)
    Xc::SVector{3} = X - p.x[1]
    return Xc
end

(cf::CalcFactor{<:PriorPoint3})(X, p::AbstractVector) = X - p

function (cf::CalcFactor{<:Point3Point3})(X, p, q)
    #
    return X[1:3] .- (q[1:3] .- p[1:3])
end
