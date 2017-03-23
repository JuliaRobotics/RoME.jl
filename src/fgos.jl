# overload for new types available via RoME

function convert{PT <: PackedInferenceType, T <: FunctorInferenceType}(::Type{PT}, ::T)
  eval(parse("$(T.name.module).Packed$(T.name.name)"))
end
function convert{T <: FunctorInferenceType, PT <: PackedInferenceType}(::Type{T}, ::PT)
  eval(parse("$(PT.name.module).$(string(PT.name.name)[7:end])"))
end
