#Work in progress

# Examples of constraint functions which can be used, however,
# \/\/\/these won't pack into ProtoBuf yet -- need to write special converters
# for Neo4j DataBase storage of complicated types. See constraints hereafter
# for more standard types.
type UniPriorPose2D <: Singleton
  Z::Distributions.MvNormal
end
type GMMPriorPose2D <: Singleton
  Z::Array{Distributions.MvNormal,1}
  W::Array{Float64,1}
end
type KDEPriorPoint2D <: Singleton
  Z::BallTreeDensity
end
type KDERangePoint2D <: Pairwise
  Z::BallTreeDensity
end
# ^^^will work on these later
