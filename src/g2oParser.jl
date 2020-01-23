# this file is for g2o integration

## Common functions for g2o parsing





## Import g2o functions






## Export g2o functions



"""
    $SIGNATURES

Export a factor graph to g2o file format.

Note:
- This funtion only supports Gaussian (i.e. Normal/MvNormal) factors, unpredictable witchcraft is used in other cases such as `AliasingScalarSampler` factor models.
"""
function exportG2o(dfg::AbstractDFG; poseRegex::Regex=r"x\d", solvable::Int=0)
  uniqVarInt = -1
  varIntLabel = Dict{Symbol, Int}()
  # all variables
  vars = ls(fg, poseRegex, solvable=solvable) |> sortDFG
  # all factors
  fcts = lsf(fg, solvable=solvable)
  # build text file based on factors, using pose variable order as guide
  for vs in vars
    # assign a unique number
    uniqVarInt += 1
    varIntLabel[vs] = uniqVarInt
    # all factors connected to variable
    vfcs = ls(fg, vs, solvable=solvable)
    kvfcs = intersect(vfcs, fcts)
    for fc in kvfcs
      isPrior(dfg, fc) ? filter!(x->x!=fc, fcts)
    end
  end

  error("Not implemented yet")
end




##
