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
function exportG2o(dfg::AbstractDFG)
  error("Not implemented yet")
end




##
