
type TopographicPoint2D <: Singleton
  Heatmap::Array{Float64,2}
  radiusOfInterest::Float64
end

function evalPotential(topo::TopographicPoint2D, Xi::Array{Graphs.ExVertex,1}; N::Int64=100)
  warn("evalPotential TopographicPoint2D is not finished yet")
  return topo.Heatmap
end
