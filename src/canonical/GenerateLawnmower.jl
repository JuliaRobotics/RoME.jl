# generate lawnmower pattern graphs



function generateCanonicalFG_Lawnmower!(vartype::IIF.InstanceType{T}
                                        numposes::Integer;
                                        dfg::AbstractDFG=initfg(),
                                        stepdistance::Real=10.0,
                                        stepspertraverse::Integer=10,
                                        stepspershift::Integer=1,
                                        turn::Symbol=:right,
                                        μ0=[0.0;0.0;pi/2] ) where {T <: InferenceVariable}
  #
  
  generateCanonicalFG_ZeroPose(; varType=T, dfg=dfg, μ0=view(μ0,1:getDimension(varType)))

  posesadded = 0
  # while posesadded < numposes
    # do traverse edge
    # change turn direction
    # do shift
    # change turn direction
  # end

  
  return dfg
end



generateCanonicalFG_Lawnmower!( w... ;kw... ) = generateCanonicalFG_Lawnmower!(Pose2, w...; kw...)
#

#