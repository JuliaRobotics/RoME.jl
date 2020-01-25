# this file is for g2o integration

export importG2o, exportG2o

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





## Export g2o functions

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc::Pose2Pose2,
                    varIntLabel::Dict,
                    uniqVarInt::Vector{Int})::String
  #
  INF = 1 ./ sqrt(fnc.z.Σ.mat)
  INF[INF .== Inf] .= 0
  varlist = Int[]
  # get variable numbers
  for vari in solverData(getFactor(dfg,fc)).fncargvID
    if !haskey(varIntLabel, vari)
      uniqVarInt[1] += 1
      varIntLabel[vari] = uniqVarInt[1]
    end
    push!(varlist, varIntLabel[vari])
  end
  return "EDGE_SE2 $(varlist[1]) $(varlist[2]) $(fnc.z.μ[1]) $(fnc.z.μ[2]) $(fnc.z.μ[3]) $(INF[1,1]) $(INF[1,1]) $(INF[1,2]) $(INF[1,3]) $(INF[2,2]) $(INF[2,3]) $(INF[2,3]) $(INF[2,3])"
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc::Pose2Point2BearingRange,
                    varIntLabel::Dict,
                    uniqVarInt::Vector{Int})::String
  #
  return "BEARRANGE not implemented"
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc,
                    varIntLabel::Dict,
                    uniqVarInt::Vector{Int})
  #
  error("unknown factor type $fnc")
end

"""
    $SIGNATURES

Export a factor graph to g2o file format.

Note:
- This funtion only supports Gaussian (i.e. Normal/MvNormal) factors, unpredictable witchcraft is used in other cases such as `AliasingScalarSampler` factor models.
"""
function exportG2o(dfg::AbstractDFG;
                   poseRegex::Regex=r"x\d",
                   solvable::Int=0,
                   ignorePriors::Bool=true,
                   filename::AbstractString="/tmp/test.txt" )
  #
  uniqVarInt = Int[-1;]
  varIntLabel = Dict{Symbol, Int}()
  # all variables
  vars = ls(dfg, poseRegex, solvable=solvable) |> sortDFG
  # all factors
  fcts = lsf(dfg, solvable=solvable)
  # build text file based on factors, using pose variable order as guide
  io = open(filename,"w")
  for vs in vars
    # assign a unique number
    # all factors connected to variable
    vfcs = ls(dfg, vs, solvable=solvable)
    kvfcs = intersect(vfcs, fcts)
    for fc in kvfcs
      # ignore priors
      ignorePriors && isPrior(dfg, fc) ? (filter!(x->x!=fc, fcts); continue) : nothing
      # only add factors to g2o file once, remove if found
      !(fc in fcts) ? continue : filter!(x->x!=fc, fcts)
      # actually add the factor to the file
      fnc = getFactorType(dfg, fc)
      # if fnc isa Pose2Pose2
        pstr = stringG2o!(dfg, fc, fnc, varIntLabel, uniqVarInt)
        println(io, pstr)
        # INF = 1 ./ sqrt(fnc.z.Σ.mat)
        # INF[INF .== Inf] .= 0
        # varlist = Int[]
        # # get variable numbers
        # for vari in solverData(getFactor(dfg,fc)).fncargvID
        #   if !haskey(varIntLabel, vari)
        #     uniqVarInt += 1
        #     varIntLabel[vari] = uniqVarInt
        #   end
        #   push!(varlist, varIntLabel[vari])
        # end
        # println(io, "EDGE_SE2 $(varlist[1]) $(varlist[2]) $(fnc.z.μ[1]) $(fnc.z.μ[2]) $(fnc.z.μ[3]) $(INF[1,1]) $(INF[1,1]) $(INF[1,2]) $(INF[1,3]) $(INF[2,2]) $(INF[2,3]) $(INF[2,3]) $(INF[2,3])" )
      # end
    end
  end
  close(io)

  return filename
end




##
