

## =============================================
## Explicit conversions
## =============================================


function homography_to_coordinates(
  ::typeof(SpecialEuclideanGroup(3; variant=:right)),
  H::AbstractMatrix{<:Real}
)
  @warn "TODO: maybe deprecate homography_to_coordinates"
  G = SpecialOrthogonalGroup(3)
  [H[1:3,4]; vee(LieAlgebra(G), log(G, H[1:3,1:3]))]
end

function coordinates_to_homography(
  ::typeof(SpecialEuclideanGroup(3; variant=:right)),
  c::AbstractVector
)
  @warn "TODO: maybe deprecate coordinates_to_homography"
  G = SpecialOrthogonalGroup(3)
  H = zeros(4, 4)
  H[1:3, 1:3] .= exp(G, hat(LieAlgebra(G), c[4:6]))
  H[1:3, 4] = c[1:3]
  H[4, 4] = 1.0
  return H
end


## =============================================
## Legacy code below
## =============================================


function convert(::Type{_Rot.QuatRotation}, q::TransformUtils.Quaternion)
  _Rot.QuatRotation(q.s, q.v...)
end
function convert(::Type{_Rot.QuatRotation}, x::SO3)
  q = convert(TransformUtils.Quaternion, x)
  convert(_Rot.QuatRotation, q)
end
function convert(::Type{T}, x::SO3) where {T <: CoordinateTransformations.AffineMap}
  LinearMap( convert(_Rot.QuatRotation, x) )
end

function convert(::Type{T}, x::SE3) where {T <: CoordinateTransformations.AffineMap}
  Translation(x.t...) ∘ convert(AffineMap{_Rot.QuatRotation{Float64}}, x.R)
end
function convert(::Type{SE3}, x::T) where {T <: CoordinateTransformations.AffineMap{_Rot.QuatRotation{Float64}}}
  SE3(x.translation[1:3], TransformUtils.Quaternion(x.linear.w, [x.linear.x,x.linear.y,x.linear.z]) )
end

function convert(::Type{SE3}, t::Tuple{Symbol, Vector{Float64}})
  if t[1]==:XYZqWXYZ
    return SE3(t[2][1:3],TransformUtils.Quaternion(t[2][4],t[2][5:7]))
  else
    error("Unknown conversion type $(t[1])")
  end
end
