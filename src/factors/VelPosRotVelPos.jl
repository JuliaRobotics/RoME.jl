

"""
$(TYPEDEF)
"""
DFG.@defObservationType VelPosRotVelPos RelativeObservation TranslationGroup(3) × TranslationGroup(3)

function (cf::CalcFactor{<:VelPosRotVelPos})(
  X_vp,
  p,  # VelPos
  q   # RotVelPos
)
  [
    X_vp.x[1] .- (q.x[2] .- p.x[1]);
    X_vp.x[2] .- (q.x[3] .- p.x[2]);
  ]
end
