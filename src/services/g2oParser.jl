# This file is for g2o integration. For more information on the actual file
# formats, refer to: https://github.com/RainerKuemmerle/g2o/wiki/File-Format

export importG2o, exportG2o, findCommand, parseG2oInstruction!

## Common functions for g2o parsing


global commands = Dict(Pose2Pose2 => :EDGE_SE2,
                       Pose2Point2BearingRange => :LANDMARK,
                       Pose3Pose3 => Symbol("EDGE_SE3:QUAT"))
#

"""
    $SIGNATURES

Hackish reverse lookup for dict[val] <= key.

Notes
- See DFG BiDictMap as better future solution (TODO)
"""
function findCommand(comm::Symbol)
  global commands
  for (key, val) in commands
    if val == comm
      return key
    end
  end
  error("$comm not found in known g2o commands, $(collect(values(commands)))")
  return nothing
end


"""
    $SIGNATURES

Read in every line from a g2o file.
"""
function importG2o(input_file::String)
    instructions = []
    # Read input file line by line.
    open(input_file) do file
        for ln in eachline(file)
            pieces = split(ln, ' ', keepempty=false)
            push!(instructions, pieces)
        end
    end
    return instructions
end

"""
    $SIGNATURES

Parses a single g2o instruction to add information to factor graph.
"""
function parseG2oInstruction!(fg::AbstractDFG,
                              instruction::Array{SubString{String},1};
                              initialize::Bool=true)
    #
    infoM6 = zeros(6,6)

    if instruction[1] == "VERTEX_SE2"
          poseId = Symbol("x", instruction[2])
          x = parse(Float64,instruction[3])
          y = parse(Float64,instruction[4])
          θ = parse(Float64,instruction[5])
          #NOTE just useing some covariance as its not included
          cov = diagm([1,1,0.1])
          addVariable!(fg, poseId, Pose2)
          if initialize
              # initVariable!(fg, poseId, MvNormal([x,y,θ],cov))
              initVariable!(fg, poseId, MvNormal([x,y,θ],cov), :parametric)
          end
    elseif instruction[1] == "EDGE_SE2"
        # Need to add a relative pose measurement between two variables.
        # Parse all of the variables starting with the symbols.
        from_pose = Symbol("x", instruction[2])
        to_pose = Symbol("x", instruction[3])

        # Parse the mean values for the relative pose measurements.
        x_state = parse(Float64, instruction[4])
        y_state = parse(Float64, instruction[5])
        theta_state = parse(Float64, instruction[6])

        # Parse the information matrix and take inverse to obtain covariance.
        info_mat = [parse(Float64, instruction[7]) parse(Float64, instruction[8])  parse(Float64, instruction[9]);
                    parse(Float64, instruction[8]) parse(Float64, instruction[10]) parse(Float64, instruction[11]);
                    parse(Float64, instruction[9]) parse(Float64, instruction[11]) parse(Float64, instruction[12])]
        cov_mat = inv(info_mat)
        # NOTE workaround to ensure cov_mat is Hermitian -- TODO use inf_mat directly for MvNormal
        cov_mat += cov_mat'
        cov_mat ./= 2

        # Make sure both variables are in the FG. Otherwise add them.
        if (from_pose in ls(fg)) == false
            addVariable!(fg, from_pose, Pose2)
        end
        if (to_pose in ls(fg)) == false
            addVariable!(fg, to_pose, Pose2)
        end

        # Add a factor between the two poses with the relative measurement.
        measurement = MvNormal([x_state; y_state; theta_state], cov_mat)
        rel_pose_factor = Pose2Pose2(measurement)
        addFactor!(fg, [from_pose, to_pose], rel_pose_factor)
    elseif instruction[1] == "EDGE_SE3:QUAT"

        # MU1 = Unitary(1,ℍ)
        # ϵU1 = identity_element(MU1)
        MSO3 = SpecialOrthogonal(3)
        ϵSO3 = identity_element(MSO3)

         # Need to add a relative pose measurement between two variables.
        # Parse all of the variables starting with the symbols.
        from_pose = Symbol("x", instruction[2])
        to_pose = Symbol("x", instruction[3])

        dt = parse.(Float64, instruction[4:6])
        qvec = parse.(Float64, instruction[[10;7:9]])
        
        dR = TU.convert(SO3, TU.Quaternion(qvec[1],qvec[2:4])).R 
        X = log(MSO3, ϵSO3, dR) 
        Xc = vee(MSO3, ϵSO3, X)

        a = parse.(Float64, instruction[11:31])
        for i=1:6
            vw = view(a, (1+(i-1)*(14-i)÷2):(i*(13-i)÷2))
            infoM6[i,i:6] = vw 
            infoM6[i:6,i] = vw
        end
        cov_mat = inv(infoM6)
        # NOTE workaround to ensure cov_mat is Hermitian -- TODO use inf_mat directly for MvNormal
        cov_mat += cov_mat'
        cov_mat ./= 2
        
        # Make sure both variables are in the FG. Otherwise add them.
        if (from_pose in ls(fg)) == false
            addVariable!(fg, from_pose, Pose3)
        end
        if (to_pose in ls(fg)) == false
            addVariable!(fg, to_pose, Pose3)
        end

        #FIXME conversion between U(1,ℍ) and SO3 covariances?
        @error "3D covariance from Quaternion is not done yet!" maxlog=1
        # dq = Manifolds.Quaternion(qvec...)

        # Add a factor between the two poses with the relative measurement.
        measurement = MvNormal([dt;Xc], cov_mat)
        rel_pose_factor = Pose3Pose3(measurement)
        addFactor!(fg, [from_pose, to_pose], rel_pose_factor)
    end
    return fg
end



## Export g2o functions

function getVariableListInts!(fct, uniqVarInt, varIntLabel)
  varlist = Int[]
  for vari in getVariableOrder(fct)
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
                    varIntLabel::OrderedDict,
                    uniqVarInt::Vector{Int};
                    overwriteMapping::Dict=Dict())::String
  #
  global commands
  # get variable numbers
  varlist = getVariableListInts!(getFactor(dfg,fc), uniqVarInt, varIntLabel)
  # get information matrix
  INF = inv(fnc.Z.Σ.mat) # 1 ./ sqrt(fnc.Z.Σ.mat) # 
  INF[INF .== Inf] .= 0
  # get command
  comm = !haskey(overwriteMapping, Pose2Pose2) ? commands[Pose2Pose2] : overwriteMapping[Pose2Pose2]
  return "$comm $(varlist[1]) $(varlist[2]) $(fnc.Z.μ[1]) $(fnc.Z.μ[2]) $(fnc.Z.μ[3]) $(INF[1,1]) $(INF[1,2]) $(INF[1,3]) $(INF[2,2]) $(INF[2,3]) $(INF[3,3])"
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc::Pose2Point2BearingRange,
                    varIntLabel::OrderedDict,
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
                    fnc::Pose3Pose3,
                    varIntLabel::OrderedDict,
                    uniqVarInt::Vector{Int};
                    overwriteMapping::Dict=Dict())::String
  #
  global commands
  # get variable numbers
  varlist = getVariableListInts!(getFactor(dfg,fc), uniqVarInt, varIntLabel)
  # get information matrix
  INF = 1 ./ sqrt(fnc.Z.Σ.mat)
  INF[INF .== Inf] .= 0

  p = DFG.getPoint(Pose3, fnc.Z.μ)
  R = p.x[2]
  Q = convert(TU.Quaternion, TU.SO3(R))

  # get command
  comm = !haskey(overwriteMapping, Pose3Pose3) ? commands[Pose3Pose3] : overwriteMapping[Pose3Pose3]
  return """$comm $(varlist[1]) $(varlist[2]) \
  $(fnc.Z.μ[1]) $(fnc.Z.μ[2]) $(fnc.Z.μ[3]) \
  $(Q.v[1]) $(Q.v[2]) $(Q.v[3]) $(Q.s) \
  $(INF[1,1]) $(INF[1,2]) $(INF[1,3]) $(INF[1,4]) $(INF[1,5]) $(INF[1,6]) \
  $(INF[2,2]) $(INF[2,3]) $(INF[2,4]) $(INF[2,5]) $(INF[2,6]) \
  $(INF[3,3]) $(INF[3,4]) $(INF[3,5]) $(INF[3,6]) \
  $(INF[4,4]) $(INF[4,5]) $(INF[4,6]) \
  $(INF[5,5]) $(INF[5,6]) \
  $(INF[6,6])"""
end

function stringG2o!(dfg::AbstractDFG,
                    fc::Symbol,
                    fnc,
                    varIntLabel::OrderedDict,
                    uniqVarInt::Vector{Int};
                    overwriteMapping::Dict=Dict{Symbol,Symbol}())
  #
  global commands
  error("unknown factor type $fnc")
end

# _writeG2oLine(
#   _, 
#   io,
#   dfg::AbstractDFG, 
#   label, 
#   i,
#   solveKey
# ) = (close(io); error("exportG2o does not support $vartype, open an issue if you would like support"))

function _writeG2oLinePose2(io, dfg::AbstractDFG, label::Symbol, i::Int, solveKey::Symbol)
  # println("trying VERTEX_SE2")
  (x,y,θ) = getPPESuggested(dfg, label, solveKey)
  write(io, "VERTEX_SE2 $i $x $y $θ\n")
end

function _writeG2oLinePose3(io, dfg::AbstractDFG, label::Symbol, i::Int, solveKey::Symbol)
  # println("WHAT IS GOING ON")
  Xc = getPPESuggested(dfg, label, solveKey)
  p = getPoint(Pose3, Xc)
  x,y,z = p.x[1]
  R = p.x[2]
  Q = convert(TU.Quaternion, TU.SO3(R))
  write(io, "VERTEX_SE3:QUAT $i $x $y $z $(Q.v[1]) $(Q.v[2]) $(Q.v[3]) $(Q.s)\n")
  return nothing
end

function _doG2oLoop(io, dfg, label, i, solveKey)
  vartype = getVariableType(dfg, label)
  typename = string(typeof(vartype).name.name)
  # FIXME, HACK, WTF https://github.com/JuliaLang/julia/issues/46871#issuecomment-1318035929
  fnc = getfield(RoME, Symbol(:_writeG2oLine, typename))
  fnc(io, dfg, label, i, solveKey)
end

function _writeG2oVertexes(
  io,
  dfg::AbstractDFG,
  varIntLabel::OrderedDict,
  solveKey::Symbol
)
  #

  # io = open(filename, "a")
  for (label,i) in pairs(varIntLabel)
    _doG2oLoop(io, dfg, label, i, solveKey)
  end
  # close(io)
  return nothing
end


function _writeG2oFactors(io, dfg, vars, fcts, varIntLabel, uniqVarInt, ignorePriors; solvable, overwriteMapping)
  for vs in vars
    # assign a unique number
    # all factors connected to variable
    vfcs = ls(dfg, vs; solvable)
    kvfcs = intersect(vfcs, fcts)
    for fc in kvfcs
      # ignore priors
      ignorePriors && isPrior(dfg, fc) ? (filter!(x->x!=fc, fcts); continue) : nothing
      # only add factors to g2o file once, remove if found
      !(fc in fcts) ? continue : filter!(x->x!=fc, fcts)
      # actually add the factor to the file
      fnc = getFactorType(dfg, fc)
      pstr = stringG2o!(dfg, fc, fnc, varIntLabel, uniqVarInt; overwriteMapping)
      println(io, pstr)
    end
  end
  # close(io)

  nothing
end

"""
    $SIGNATURES

Export a factor graph to g2o file format.

Note:
- This funtion only supports Gaussian (i.e. Normal/MvNormal) factors, unpredictable witchcraft is used in other cases such as `AliasingScalarSampler` factor models.
"""
function exportG2o(
  dfg::AbstractDFG;
  poseRegex::Regex=r"x\d",
  solvable::Int=0,
  ignorePriors::Bool=true,
  filename::AbstractString="/tmp/test.txt",
  overwriteMapping::Dict=Dict{Symbol, Symbol}(),
  varIntLabel::OrderedDict{Symbol, Int} = OrderedDict{Symbol, Int}(),
  solveKey::Union{Nothing,Symbol}=nothing,
)
  #
  uniqVarInt = Int[-1;]
  # all variables
  vars = ls(dfg, poseRegex; solvable) |> sortDFG
  # all factors
  fcts = lsf(dfg; solvable)
  # build text file based on factors, using pose variable order as guide
  io = open(filename,"w")
  
  if !isnothing(solveKey)
    _writeG2oVertexes( io,dfg, varIntLabel,solveKey)
  end
  _writeG2oFactors(io, dfg, vars, fcts, varIntLabel, uniqVarInt, ignorePriors; overwriteMapping, solvable)
  close(io)

  return filename
end




##
