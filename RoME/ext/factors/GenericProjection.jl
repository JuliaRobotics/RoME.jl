

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
  fromtype::String
  totype::String
  cam::Dict{Symbol,Any}
  Z::PackedSamplableBelief
end
function convert(::Type{GenericProjection}, packed::PackedGenericProjection)
  fromtype = DFG.getTypeFromSerializationModule(packed.fromtype) # Pose3
  totype   = DFG.getTypeFromSerializationModule(packed.totype)   # Point3
  cam = convert(DFG.getTypeFromSerializationModule(packed.cam[:_type]), packed.cam)
  Z = convert(SamplableBelief, packed.Z) 
  return GenericProjection{fromtype,totype}(
      cam,
      Z
    )
  end
  function convert(::Type{PackedGenericProjection}, obj::GenericProjection{FT,TT}) where {FT,TT}
    camdict = Dict{Symbol, Any}(
      :height => obj.cam.height,
      :width => obj.cam.width,
      :kc => obj.cam.kc,
      :K => obj.cam.K[:],
      :Ki => obj.cam.Ki[:],
      :_type => "CameraModels.CameraCalibration"
    )
  return PackedGenericProjection(
    string(FT.name.module, ".", FT.name.name), # "RoME.Pose3",  # string(FT)
    string(TT.name.module, ".", TT.name.name), # "RoME.Point3", # string(TT)
    camdict,
    convert(PackedSamplableBelief, obj.Z) 
  )
end