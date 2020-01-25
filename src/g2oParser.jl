# this file is for g2o integration

export importG2o, exportG2o

## Common functions for g2o parsing


global commands = Dict(Pose2Pose2 => :EDGE_SE2,
                       Pose2Point2BearingRange => :LANDMARK)


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

function getVariableListInts!(fct, uniqVarInt, varIntLabel)
  varlist = Int[]
  for vari in solverData(fct).fncargvID
    if !haskey(varIntLabel, vari)
      uniqVarInt[1] += 1
      varIntLabel[vari] = uniqVarInt[1]
    end
    push!(varlist, varIntLabel[vari])
  end
  varlist
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc::Pose2Pose2,
                    varIntLabel::Dict,
                    uniqVarInt::Vector{Int};
                    overwriteMapping::Dict=Dict())::String
  #
  global commands
  # get variable numbers
  varlist = getVariableListInts!(getFactor(dfg,fc), uniqVarInt, varIntLabel)
  # get information matrix
  INF = 1 ./ sqrt(fnc.z.Σ.mat)
  INF[INF .== Inf] .= 0
  # get command
  comm = !haskey(overwriteMapping, Pose2Pose2) ? commands[Pose2Pose2] : overwriteMapping[Pose2Pose2]
  return "$comm $(varlist[1]) $(varlist[2]) $(fnc.z.μ[1]) $(fnc.z.μ[2]) $(fnc.z.μ[3]) $(INF[1,1]) $(INF[1,1]) $(INF[1,2]) $(INF[1,3]) $(INF[2,2]) $(INF[2,3]) $(INF[2,3]) $(INF[2,3])"
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc::Pose2Point2BearingRange,
                    varIntLabel::Dict,
                    uniqVarInt::Vector{Int};
                    overwriteMapping::Dict=Dict{Symbol,Symbol}())::String
  #
  global commands
  # get variable numbers
  varlist = getVariableListInts!(getFactor(dfg,fc), uniqVarInt, varIntLabel)
  # get information matrix
  INF = [1/sqrt(fnc.bearing.σ); 0; 1/sqrt(fnc.range.σ)]
  # get command
  comm = !haskey(overwriteMapping, Pose2Point2BearingRange) ? commands[Pose2Point2BearingRange] : overwriteMapping[Pose2Point2BearingRange]
  return "$comm $(varlist[1]) $(varlist[2]) $(fnc.bearing.μ) $(fnc.range.μ) $(INF[1]) $(INF[2]) $(INF[3])"
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc,
                    varIntLabel::Dict,
                    uniqVarInt::Vector{Int};
                    overwriteMapping::Dict=Dict{Symbol,Symbol}())
  #
  global commands
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
                   filename::AbstractString="/tmp/test.txt",
                   overwriteMapping::Dict=Dict{Symbol, Symbol}())
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
      pstr = stringG2o!(dfg, fc, fnc, varIntLabel, uniqVarInt, overwriteMapping=overwriteMapping)
      println(io, pstr)
    end
  end
  close(io)

  return filename
end




##
