
"""
$(TYPEDEF)

Introduce direct observations on all dimensions of a Pose2 variable:

Example:
--------
```julia
PriorVelPos3( MvNormal(zeros(6), Matrix(Diagonal(ones(6).^2))) )
```
"""
DFG.@defObservationType PriorVelPos3 PriorObservation TranslationGroup(3) ×
                                                      TranslationGroup(3)

function (cf::CalcFactor{<:PriorVelPos3})(m::ArrayPartition, p::ArrayPartition)
    M = getManifold(PriorVelPos3)
    # TODO, Lie Group for now, expand to Riemannian
    ε = getPointIdentity(M)
    Xc = vee(M, ε, log(M, p, m))
    return Xc
end

#TODO Serialization of reference point p 
## Serialization support
