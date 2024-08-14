

"""
$(TYPEDEF)
"""
Base.@kwdef struct VelAlign{D <: IIF.SamplableBelief} <: AbstractManifoldMinimize
  """ Measurement vector of [1;0;0] here implies body velocity is assumed forward along body x-axis """
  Z::D = MvNormal([1;0;0.0],diagm([0.1;0.1;0.1]))
end

getManifold(::InstanceType{VelAlign}) = TranslationGroup(3)



function IncrementalInference.preambleCache(
  dfg::AbstractDFG, 
  vars::AbstractVector{<:DFGVariable}, 
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
  p_V = norm(cf.cache.p_vel(w_T_p)) .* X_v
  q_V = cf.cache.q_rot(w_T_q)' * cf.cache.q_vel(w_T_q)
  p_V - p_R_q * q_V
end


"""
$(TYPEDEF)

Serialization type for `VelAlign`.
"""
Base.@kwdef struct PackedVelAlign <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{VelAlign}, d::PackedVelAlign)
  return VelAlign( convert(SamplableBelief, d.Z) )
end
function convert(::Type{PackedVelAlign}, d::VelAlign)
  return PackedVelAlign( convert(PackedSamplableBelief, d.Z) )
end