

"""
$(TYPEDEF)
"""
Base.@kwdef struct VelPosRotVelPos{D <: IIF.SamplableBelief} <: AbstractManifoldMinimize
    Z::D = MvNormal(zeros(3),diagm([0.1;0.1;0.1]))
end

getManifold(::InstanceType{VelPosRotVelPos}) = ProductGroup(
  ProductManifold(
    TranslationGroup(3),
    TranslationGroup(3)
  )
)

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


"""
$(TYPEDEF)

Serialization type for `VelPosRotVelPos`.
"""
Base.@kwdef struct PackedVelPosRotVelPos <: AbstractPackedFactor
    Z::PackedSamplableBelief
end
function convert(::Type{VelPosRotVelPos}, d::PackedVelPosRotVelPos)
  return VelPosRotVelPos( convert(SamplableBelief, d.Z) )
end
function convert(::Type{PackedVelPosRotVelPos}, d::VelPosRotVelPos)
  return PackedVelPosRotVelPos( convert(PackedSamplableBelief, d.Z) )
end