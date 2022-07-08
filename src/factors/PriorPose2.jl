
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
  Z::T = MvNormal(zeros(3), diagm([1;1;0.1]))
end

DFG.getManifold(::InstanceType{PriorPose2}) = getManifold(Pose2) # SpecialEuclidean(2)

function (cf::CalcFactor{<:PriorPose2})(m, p)
  M = getManifold(Pose2)
  Xc = vee(M, p, log(M, p, m))
  return Xc

  # M = getManifold(Pose2)
  # # X = allocate(p)
  # # X = ProductRepr(zeros(MVector{2}), zeros(MMatrix{2,2}))
  # # log!(M, X, p, m)
  # X = log(M, p, m)
  # return  SA[X.x[1][1],X.x[1][2],X.x[2][2]]

end
# BenchmarkTools.Trial: 10000 samples with 1 evaluation.
#  Range (min … max):  22.299 μs …  11.920 ms  ┊ GC (min … max): 0.00% … 98.28%
#  Time  (median):     29.930 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   40.097 μs ± 171.479 μs  ┊ GC (mean ± σ):  5.60% ±  1.39%

#    ▂▆█▇▆▅▄▃▃▄▃▂        ▂▂▂▂▁                                   ▂
#   ▆██████████████▆▆▆▄▅▇██████▇▇█▆▇▆▄▅▄▆▆▅▅▅▆▅▅▅▄▅▄▄▄▅▅▅▄▄▄▃▅▄▄ █
#   22.3 μs       Histogram: log(frequency) by time       127 μs <

#  Memory estimate: 26.88 KiB, allocs estimate: 481.

#TODO Serialization of reference point p 
## Serialization support

"""
$(TYPEDEF)
"""
Base.@kwdef struct PackedPriorPose2  <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{PackedPriorPose2}, d::PriorPose2)
  return PackedPriorPose2(convert(PackedSamplableBelief, d.Z))
end
function convert(::Type{PriorPose2}, d::PackedPriorPose2)
  return PriorPose2(convert(SamplableBelief, d.Z))
end




## NOTE likely deprecated comparitors, see DFG compareFields, compareAll instead
function compare(a::PriorPose2,b::PriorPose2; tol::Float64=1e-10)
  compareDensity(a.Z, b.Z)
end
