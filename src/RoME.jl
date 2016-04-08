module RoME

export 
  measureMeanDist,
  predictBodyBR,
  getLastPose2D,
  odomKDE,
  addOdoFG!,
  newLandm!,
  addBRFG!,
  addMMBRFG!,
  projNewLandm!,
  malahanobisBR,
  initFactorGraph!,
  
  # Some vizualization tools
  plotLsrScanFeats,
  drawFeatTrackers,
  saveImgSeq,
  lsrBR,
  
  

include("robots/BayesTracker.jl")
include("robots/RobotViz.jl")
include("robots/RobotUtils.jl")
include("robots/SimulationUtils.jl")
  
end


