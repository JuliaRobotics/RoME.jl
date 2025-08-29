
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorPose2( MvNormal([10; 10; pi/6.0], Matrix(Diagonal([0.1;0.1;0.05].^2))) )
```
"""
Base.@kwdef struct PriorPose2{T <: SamplableBelief} <: IIF.AbstractPrior
    Z::T = MvNormal(zeros(3), diagm([1; 1; 0.1]))
end

DFG.getManifold(::InstanceType{PriorPose2}) = getManifold(Pose2) # SpecialEuclidean(2)

function (cf::CalcFactor{<:PriorPose2})(
    _m::AbstractArray{MT},
    _p::AbstractArray{PT},
) where {MT <: Real, PT <: Real}
    T = promote_type(MT, PT)
    m = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _m)
    p = convert(ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}}, _p)
    return cf(m, p)
end

function (cf::CalcFactor{<:PriorPose2})(
    m::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
    p::ArrayPartition{T, Tuple{SVector{2, T}, SMatrix{2, 2, T, 4}}},
) where {T <: Real}
    M = getManifold(PriorPose2)
    X = log(M, p, m) # Currently X ∈ TₚM, #TODO should it be TₘM? Also update the rest if this is wrong.
    return vee(LieAlgebra(M), X)
end

#TODO Serialization of reference point p 
## Serialization support

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPriorPose2 <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function DFG.pack(d::PriorPose2)
    return PackedPriorPose2(packDistribution(d.Z))
end
function DFG.unpack(d::PackedPriorPose2)
    return PriorPose2(unpackDistribution(d.Z))
end

## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
function compare(a::PriorPose2, b::PriorPose2; tol::Float64 = 1e-10)
    return compareDensity(a.Z, b.Z)
end
