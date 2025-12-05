# 2D SLAM with velocity states

"""
$(TYPEDEF)
"""
DFG.@defObservationType DynPoint2VelocityPrior PriorObservation TranslationGroup(4)

"""
$(TYPEDEF)
"""
DFG.@defObservationType DynPoint2DynPoint2 RelativeObservation TranslationGroup(4)

function (cfo::CalcFactor{<:DynPoint2DynPoint2})(z, xi, xj)
    #
    dt = DFG.calcDeltatime(cfo.fullvariables[1], cfo.fullvariables[2])
    res12 = z[1:2] - (xj[1:2] - (xi[1:2] + dt * xi[3:4]))
    res34 = z[3:4] - (xj[3:4] - xi[3:4])
    return [res12; res34]
end

"""
$(TYPEDEF)
"""
DFG.@defObservationType Point2Point2Velocity RelativeObservation TranslationGroup(4)

function (cfo::CalcFactor{<:Point2Point2Velocity})(z, xi, xj)
    #
    dt = DFG.calcDeltatime(cfo.fullvariables[1], cfo.fullvariables[2])
    dp = (xj[1:2] .- xi[1:2])
    dv = (xj[3:4] .- xi[3:4])

    res12 = z[1:2] .- dp
    res34 = dp / dt .- 0.5 * (xj[3:4] .+ xi[3:4])  # (dp/dt - 0.5*(xj[3:4]+xi[3:4])) # midpoint integration

    return [res12; res34]
end

#
