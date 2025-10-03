
struct MixtureFluxPose2Pose2{F <: AbstractObservation} <: AbstractRelativeMinimize
    Z::F
    # delta time between variables
    DT::Base.RefValue{Float64}
end

mutable struct PackedMixtureFluxPose2Pose2 <: AbstractPackedObservation
    Z::PackedMixture
    DT::Float64
end
