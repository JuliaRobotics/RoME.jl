

"""
$(TYPEDEF)
"""
DFG.@defObservationType VelAlign RelativeObservation TranslationGroup(3)

function IncrementalInference.preambleCache(
  dfg::AbstractDFG, 
  vars::AbstractVector{<:VariableCompute}, 
  ::VelAlign,
)
  # TODO, obsolete -- replace with NamedTuple submanifold checks
  @assert typeof(getVariableType(vars[1])) <: VelPos3 "VelAlign expects first variable type VelPos3"
  @assert typeof(getVariableType(vars[2])) <: RotVelPos "VelAlign expects second variable type RotVelPos"
  @assert typeof(getVariableType(vars[3])) <: Rotation3 "VelAlign expects third variable type Rotation"
  (;
    p_vel= s->s.x[1],
    q_rot= s->s.x[1],
    q_vel= s->s.x[2],
  )
end

function (cf::CalcFactor{<:VelAlign})(
  X_v,
  w_T_p,  # VelPos
  w_T_q,  # RotVelPos
  p_R_q
)
  # body velocity scaled by real speed
  # FIXME align directions and adaptively weight by speed (i.e. slow is less important)
  p_V = norm(cf.cache.p_vel(w_T_p)) .* X_v
  q_V = cf.cache.q_rot(w_T_q)' * cf.cache.q_vel(w_T_q)
  p_V - p_R_q * q_V
end
