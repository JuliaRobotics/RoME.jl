

"""
    $TYPEDEF

Generic projections, e.g. camera making sparse feature sighting:

```julia
camcal = CameraModels.PinholeCamera(480,640,410,410,241,319) # simple model example
sigma = [10 0; 0 10.0]
meas = (pixel_row,pixel_col)
gp = GenericProjection{
  Pose3,
  Point3
} = (
  MvNormal(
    meas,
    sigma
  ),
  camcal
)
```
"""
@kwdef struct GenericProjection{SRC<:InferenceVariable,TRG<:InferenceVariable,C,D} <: AbstractManifoldMinimize
  cam::C
  Z::D
end

GenericProjection{SRC,TRG}(cam::C, Z::D) where {SRC<:InferenceVariable,TRG<:InferenceVariable,C,D} = GenericProjection{SRC,TRG,C,D}(;cam,Z)


getManifold(gp::GenericProjection{S,T}) where {S, T} = TranslationGroup(getDimension(gp.Z))



"""
$(TYPEDEF)

Serialization type for `GenericProjection`.
"""
Base.@kwdef struct PackedGenericProjection <: AbstractPackedFactor
  cam::Dict{Symbol,Any}
  Z::PackedSamplableBelief
end
function convert(::Type{GenericProjection}, packed::PackedGenericProjection)
  cam = convert(DFG.getTypeFromSerializationModule(packed.cam[:_type]), packed.cam)
  Z = convert(SamplableBelief, packed.Z) 
  # FIXME, dont hardcode Pose3,Point3
  return GenericProjection{Pose3,Point3}( 
    cam,
    Z
  )
end
function convert(::Type{PackedGenericProjection}, obj::GenericProjection)
  camdict = Dict{Symbol, Any}(
    :height => obj.cam.height,
    :width => obj.cam.width,
    :kc => obj.cam.kc,
    :K => obj.cam.K[:],
    :Ki => obj.cam.Ki[:],
    :_type => "CameraModels.CameraCalibration"
  )
  return PackedGenericProjection( 
    camdict,
    convert(PackedSamplableBelief, obj.Z) 
  )
end