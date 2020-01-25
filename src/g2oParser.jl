# This file is for g2o integration. For more information on the actual file
# formats, refer to: https://github.com/RainerKuemmerle/g2o/wiki/File-Format

## Common functions for g2o parsing





## Import g2o functions
function importG2o(input_file::String)
    instructions = []
    # Read input file line by line.
    open(input_file) do file
        for ln in eachline(file)
            println("$(length(ln)), $(ln)")
            pieces = split(ln, ' ', keepempty=false)
            @show pieces
            push!(instructions, pieces)
        end
    end
    return instructions
end

function parseG2oInstruction(instruction::Array{SubString{String},1})
    #
end



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
