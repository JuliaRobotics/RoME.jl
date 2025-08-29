

struct MixtureFluxPose2Pose2{F <: AbstractFactor} <: AbstractRelativeMinimize
  Z::F
  # delta time between variables
  DT::Base.RefValue{Float64}
end

mutable struct PackedMixtureFluxPose2Pose2 <: AbstractPackedFactor
  Z::PackedMixture
  DT::Float64
end