

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
@kwdef struct GenericProjection{SRC<:InferenceVariable,TRG<:InferenceVariable,C,D}
  Z::D
  calib::C
end

getManifold(gp::GenericProjection{S,T}) where {S, T} = TranslationGroup(getDimension(gp.Z))
