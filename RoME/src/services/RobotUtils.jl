

"""
    $SIGNATURES

Remove all marginalization and make solvable (=1) all variables and factors that don't contain `drt` in their label.

Notes
- Dead Reckon Tether (DRT)

DevNotes
- Legacy, poor implementation, needs to be improved

Related

dontMarginalizeVariablesAll!, setSolvable!, defaultFixedLagOnTree!
"""
function enableSolveAllNotDRT!(dfg::AbstractDFG; solvable::Int=1)
  dontMarginalizeVariablesAll!(dfg)
  foreach(x->setSolvable!(dfg, x, solvable), ls(dfg))
  foreach(x->setSolvable!(dfg, x, solvable), lsf(dfg))
  nothing
end


"""
    $(SIGNATURES)

Calculate the cartesian distance between two vertices in the graph using their symbol name, and by maximum belief point.
"""
function getRangeKDEMax2D(fgl::AbstractDFG, vsym1::Symbol, vsym2::Symbol)
  x1 = getKDEMax(getBelief(fgl, vsym1))
  x2 = getKDEMax(getBelief(fgl, vsym2))
  norm(x1[1:2]-x2[1:2])
end


"""
    $SIGNATURES

Return the last `number::Int` of poses according to `filterLabel::Regex`.

Notes
- Uses FIFO add history of variables to the distribued factor graph object as search index.

DevNotes
- TODO consolidate with generic filter by regext on variable labels in IIF instead.
"""
function getLastPoses(dfg::AbstractDFG;
                      filterLabel::Regex=r"x\d",
                      number::Int=5)::Vector{Symbol}
  #
  # filter according to pose label
  syms = filter(l->occursin(filterLabel, string(l)), getAddHistory(dfg))

  # return the last segment of syms
  len = length(syms)
  st = number < len ? len-number+1 : 1
  return syms[st:end]
end

"""
    $SIGNATURES

Set old poses and adjacent factors to `solvable::Int=0` (default).

Notes
- `youngest::Int` and `oldest::Int` set the limits of search by count,
  - `oldest` set large enough for solver loop defintely disengage old parts in re-occuring cycle.
- `filterLabel::Regex` sets the template to search for each pose label.
- Poses are assumed to be a thread through time that connects the local exploration variables.
- Initially developed to remove old variables and factors from a solution, in combination with fix-lag marginalization.
  - `getSolverParams(fg).isfixedlag=true`

Related:

getLastPoses
"""
function setSolvableOldPoses!(dfg::AbstractDFG;
                              youngest::Int=50,
                              oldest::Int=200,
                              solvable=0,
                              filterLabel::Regex=r"x\d")
  #
  # collect old variables and factors to disable from next solve
  newPoses = getLastPoses(dfg,filterLabel=filterLabel, number=youngest)
  oldPoses = setdiff(getLastPoses(dfg,filterLabel=filterLabel, number=oldest), newPoses)
  allFcts = (oldPoses .|> x->ls(dfg,x))
  fctAdj = 0 < length(allFcts) ? union(allFcts...) : Symbol[]

  # all together
  fullList = [oldPoses; fctAdj]

  # use solvable=0 to disable variables and factors in the next solve
  map(x->setSolvable!(dfg, x, solvable), fullList)

  return fullList
end



"""
    $(SIGNATURES)

Initialize a factor graph object as Pose2, Pose3, or neither and returns variable and factor symbols as array.
"""
function initFactorGraph!(fg::AbstractDFG;
                          P0::Union{Array{Float64,2},Nothing}=nothing,
                          init::Union{Vector{Float64},Nothing}=nothing,
                          N::Int=100,
                          lbl::Symbol=:x0,
                          solvable::Int=1,
                          firstPoseType=Pose2,
                          labels::Vector{Symbol}=Symbol[])
  #
  nodesymbols = Symbol[]
  if firstPoseType == Pose2
      init = init!=nothing ? init : zeros(3)
      P0 = P0!=nothing ? P0 : Matrix(Diagonal([0.03;0.03;0.001]))
      # init = vectoarr2(init)
      addVariable!(fg,lbl,Pose2,N=N, solvable=solvable, tags=labels )
      push!(nodesymbols, lbl)
      # v1 = addVariable!(fg, lbl, init, P0, N=N, solvable=solvable, tags=labels)
      fctVert = addFactor!(fg, [lbl;], PriorPose2(MvNormal(init, P0)), solvable=solvable, tags=labels) #[v1],
      push!(nodesymbols, Symbol(fctVert.label))
  end
  if firstPoseType == Pose3
      init = init!=nothing ? init : zeros(6)
      P0 = P0!=nothing ? P0 : Matrix(Diagonal([0.03;0.03;0.03;0.001;0.001;0.001]))
      addVariable!(fg,lbl,Pose2,N=N,solvable=solvable,tags=labels )
      push!(nodesymbols, lbl)
      # v1 = addVariable!(fg, lbl, init, P0, N=N, solvable=solvable, tags=labels)
      fctVert = addFactor!(fg, [lbl;], PriorPose3(MvNormal(init, P0)), solvable=solvable, tags=labels) #[v1],
      push!(nodesymbols, Symbol(fctVert.label))
  end
  return nodesymbols
end



function replaceFactorPose3Pose3Mean!(
  dfg::AbstractDFG,
  flb::Symbol,
  H::AbstractMatrix;
  graphinit=false
)
  fct = getFactor(dfg, flb)
  vars = getVariableOrder(fct)
  mn, sig = getMeasurementParametric(getFactorType(fct))
  mn_ = homography_to_coordinates(getManifold(Pose3), float.(H))
  @info "Replace factor with new mean" string(mn') string(mn_')
  tags = IIF.getTags(fct) |> collect
  deleteFactor!(dfg, flb)
  addFactor!(
    dfg, 
    vars, 
    Pose3Pose3(
      MvNormal(
        mn_,
        sig
      )
    );
    tags,
    graphinit
  )
end



# ------------------------------------
# Transfered from IncrementalInference



function get2DSamples(fg::AbstractDFG;
                      from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),
                      minnei::Int=0,
                      regexKey::Regex=r"x")
  #
  X = Vector{Float64}()
  Y = Vector{Float64}()

  # if sym = 'l', ignore single measurement landmarks
  allids = listVariables(fg, regexKey)  # fg.IDs
  saids = DFG.sortDFG(allids)
  for id in saids
    vertlbl = string(id)
    val = parse(Int,split(vertlbl[2:end],'_')[1])
    if from <= val && val <= to
      if length( listNeighbors(fg, id ) ) >= minnei
        # if length(out_neighbors(fg.v[id[2]],fg.g)) >= minnei
        M = getManifold(getVariable(fg, id))
        pts = getPoints(fg, id)
        pts_ = AMP.makeCoordsFromPoint.(Ref(M), pts)
        @cast X_[i] := pts_[i][1]
        @cast Y_[i] := pts_[i][2]
        X = [X; X_]
        Y = [Y; Y_]
      end
    end
  end

  return X,Y
end

"""
    $SIGNATURES

List all variables that fall in numerical range `from`, `to`, and with prefix key as specified.

Related

DFG.getVariableLabelNumber, DFT.findFactorsBetweenNaive
"""
function listVariablesLabelsWithinRange(fg::AbstractDFG,
                                        regexKey::Regex=r"x";
                                        from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),
                                        minnei::Int=0)
  #

  # if sym = 'l', ignore single measurement landmarks
  allids = listVariables(fg, regexKey)  # fg.IDs
  saids = DFG.sortDFG(allids)
  mask = Array{Bool,1}(undef, length(saids))
  fill!(mask, false)
  count = 0
  for id in saids
    count += 1
    if length( listNeighbors(fg, id) ) >= minnei
      mask[count] = true
    end
    if occursin(regexKey, string(id)) && (from != 0 || to != (2^(Sys.WORD_SIZE-1)-1))
      vertlbl = string(id)
        # TODO won't work with nested labels
        m = match(r"\d+", string(id))
        val = parse(Int, m.match)
        # val_ = split(vertlbl[2:end],'_')[1]
        # val = parse(Int,val_)
        if !(from <= val <= to)
          mask[count] = false
        end
    end
  end

  saids[mask]
end


function get2DSampleMeans(fg::AbstractDFG,
                          regexKey::Regex=r"x";
                          from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),
                          minnei::Int=0)
  #
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]

  vsyms = listVariablesLabelsWithinRange(fg, regexKey, from=from, to=to, minnei=minnei)

  for id in vsyms
    X=[X; Statistics.mean( vec( getVal(fg, id )[1,:] ) )]
    Y=[Y; Statistics.mean( vec( getVal(fg, id )[2,:] ) )]
    # crude test for pose TODO probably not going to always work right
    if string(id)[1] == 'x'
      Th=[Th; Statistics.mean( vec( getVal(fg, id )[3,:] ) )]
    end
    push!(LB, string(id))
  end
  return X,Y,Th,LB
end

#draw landmark positions
function getAll2DMeans(fg, sym::Regex)
  return get2DSampleMeans(fg, sym )
end

function getAll2DPoses(fg::AbstractDFG; regexKey=r"x")
    return getAll2DSamples(fg, regexKey=regexKey )
end

function get2DPoseSamples(fg::AbstractDFG; from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1), regexKey=r"x")
  return get2DSamples(fg, regexKey=regexKey, from=from, to=to )
end

function get2DPoseMeans(fg::AbstractDFG; from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1), regexKey=r"x")
  return get2DSampleMeans(fg, regexKey, from=from, to=to )
end


function get2DPoseMax(fgl::G;
					  regexKey::Regex=r"x",
                      from::Int=-(2^(Sys.WORD_SIZE-1)-1), to::Int=(2^(Sys.WORD_SIZE-1)-1) ) where G <: AbstractDFG
  #
  # xLB,ll = ls(fgl) # TODO add: from, to, special option 'x'
  xLB = DFG.getVariableIds(fgl, regexKey)
  saids = DFG.sortDFG(xLB)
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]
  for slbl in saids
    lbl = string(slbl)
    if from <= parse(Int,split(lbl[2:end],'_')[1]) <=to
      mv = getKDEMax(getBelief(fgl,slbl))
      push!(X,mv[1])
      push!(Y,mv[2])
      push!(Th,mv[3])
      push!(LB, string(lbl))
    end
  end
  return X, Y, Th, LB
end


function get2DLandmSamples(fg::G;
                           from::Int=0,
                           to::Int=(2^(Sys.WORD_SIZE-1)-1),
                           minnei::Int=0 ) where G <: AbstractDFG
  #
  return get2DSamples(fg, regexKey=r"l", from=from, to=to, minnei=minnei )
end

function get2DLandmMeans(fg::G;
                         from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),
                         minnei::Int=0,
                         landmarkRegex::Regex=r"l" ) where G <: AbstractDFG
  #
  return get2DSampleMeans(fg, landmarkRegex, from=from, to=to, minnei=minnei )
end

function removeKeysFromArr(fgl::G,
                           torm::Array{Int,1},
                           lbl::Array{String,1}) where G <: AbstractDFG
  #
  retlbs = String[]
  for i in 1:length(lbl)
    id = parse(Int,split(lbl[i][2:end],'_')[1])
    if something(findfirst(isequal(id), torm), 0) == 0 #findfirst(torm,id) == 0
      push!(retlbs, lbl[i])
    else
      println("removeKeysFromArr -- skipping $(lbl[i]), id=$(id)")
    end
  end
  return retlbs
end
function removeKeysFromArr(fgl::G,
                           torm::Array{Int,1},
                           lbl::Array{Symbol,1} ) where G
  #
  removeKeysFromArr(fgl, torm, string.(lbl))
end

function get2DLandmMax(fgl::G;
                       from::Int=-(2^(Sys.WORD_SIZE-1)-1),
                       to::Int=(2^(Sys.WORD_SIZE-1)-1),
                       showmm=false, MM::Dict{Int,T}=Dict{Int,Int}(),
                       regexLandmark::Regex=r"l"  ) where {G <: AbstractDFG, T}
  #
  # xLB,lLB = ls(fgl) # TODO add: from, to, special option 'x'
  lLB = DFG.getVariableIds(fgl, regexLandmark)
  if !showmm lLB = removeKeysFromArr(fgl, collect(keys(MM)), lLB); end
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]
  for lb in lLB
    @show lb
    lbl = string(lb)
    if from <= parse(Int,split(lbl[2:end],'_')[1]) <=to
      mv = getKDEMax(getBelief(fgl, Symbol(lb)))
      push!(X,mv[1])
      push!(Y,mv[2])
      push!(LB, string(lbl))
    end
  end
  return X, Y, Th, LB
end



# convenience function to add DIDSON sonar constraints to graph
function addLinearArrayConstraint(fgl::G,
                                  rangebearing::Union{Tuple{Float64, Float64}, Vector{Float64}},
                                  pose::Symbol,
                                  landm::Symbol ;
                                  rangecov::Float64=3e-4,
                                  bearingcov::Float64=3e-4 ) where G <: AbstractDFG

  #

  cl = LinearRangeBearingElevation((rangebearing[1],rangecov),(rangebearing[2],bearingcov))
  if !haskey(fgl.IDs, landm)
    pts = getVal(fgl, pose) + cl
    N = size(pts,2)
    vl1 = addVariable!(fgl, landm,  pts,  N=N)
    println("Automatically added $(landm) to the factor graph")
  end
  addFactor!(fgl, [getVert(fgl, pose); getVert(fgl, landm)], cl)
  nothing
end


function addSoftEqualityPoint2D(fgl::G,
                                l1::Symbol,
                                l2::Symbol;
                                dist=MvNormal([0.0;0.0],Matrix{Float64}(LinearAlgebra.I, 2,2)),
                                solvable::Int=1  )  where G <: AbstractDFG
  #
  pp = Point2DPoint2D(dist)
  addFactor!(fgl, [l1,l2], pp, solvable=solvable)
end
