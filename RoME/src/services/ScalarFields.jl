# ScalarField related functions loaded when Interpolations.jl is available.


"""
    $SIGNATURES

Load gridded elevation data `dem` into a factor graph `fg` as a collection
of Point3 variables. Each variable is connected to its 4-neighborhood by
relative Point3Point3 MvNormal constraints with mean defined by their
relative position on the grid (`x`, `y`) and covariance `meshEdgeSigma`.
"""
function _buildGraphScalarField!( fg::AbstractDFG,
                                  dem::AbstractMatrix, # assume grayscale image for now
                                  x::AbstractVector,
                                  y::AbstractVector;
                                  solvable::Int=0,
                                  marginalized::Bool=true,
                                  meshEdgeSigma=diagm([1;1;1]),
                                  refKey::Symbol = :simulated )
  #
  # assume regular grid
  dx, dy = x[2]-x[1], y[2]-y[1]
  for i in 1:length(x)
    for j in 1:length(y)
      s = Symbol("pt$(i)_$(j)") # unique identifier
      pt = addVariable!(fg, s, Point3, solvable=solvable)    
      setMarginalized!(pt, marginalized) # assume solveKey=:default
      
      # ...
      refVal = [x[i];y[j];dem[i,j]]
      simPPE = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
      setPPE!(pt, refKey, typeof(simPPE), simPPE)
      
      # Regular grid triangulation:
      #  add factor to (i-1,j)   |
      #  add factor to (i, j-1)  -
      #  add factor to (i-1, j-1)  \

      # no edges to prev row on first row
      if i>1
        dVal1 = dem[i,j]-dem[i-1,j]
        f = Point3Point3(MvNormal([dx, 0, dVal1], meshEdgeSigma))
        addFactor!(fg, [Symbol("pt$(i-1)_$(j)"), s], f, solvable=solvable, graphinit=false)
      end
      
      # no edges to prev column on first column
      if j>1
        dVal2 = dem[i,j]-dem[i,j-1]
        f = Point3Point3(MvNormal([0, dy, dVal2], meshEdgeSigma))
        addFactor!(fg, [Symbol("pt$(i)_$(j-1)"),s], f, solvable=solvable, graphinit=false)
      end

      # no edges to add on first element
      if i>1 && j>1
        dVal3 = dem[i,j]-dem[i-1,j-1]
        f = Point3Point3(MvNormal([dx, dy, dVal3], meshEdgeSigma))
        addFactor!(fg,[Symbol("pt$(i-1)_$(j-1)"),s], f, solvable=solvable, graphinit=false)
      end
      
    end
  end

  nothing
end




#