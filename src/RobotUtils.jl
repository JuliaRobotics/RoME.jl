
export getVariablesLabelsWithinRange, enableSolveAllNotDRT!

"""
    $SIGNATURES

Remove all marginalization and make solvable (=1) all variables and factors that don't contain `drt` in their label.

Notes
- Dead Reckon Tether (DRT)

Related

dontMarginalizeVariablesAll!, setSolvable!, defaultFixedLagOnTree!
"""
function enableSolveAllNotDRT!(dfg::AbstractDFG; solvable::Int=1)
  dontMarginalizeVariablesAll!(dfg)
  foreach(x->setSolvable!(dfg, x, solvable), ls(dfg))
  foreach(x->setSolvable!(dfg, x, solvable), lsf(dfg))
  nothing
end

mutable struct RangeAzimuthElevation
  range::Float64
  azimuth::Float64
  elevation::Union{Nothing,Float64}
end

function convert(::Type{_Rot.QuatRotation}, q::TransformUtils.Quaternion)
  _Rot.QuatRotation(q.s, q.v...)
end
function convert(::Type{_Rot.QuatRotation}, x::SO3)
  q = convert(TransformUtils.Quaternion, x)
  convert(_Rot.QuatRotation, q)
end
function convert(::Type{T}, x::SO3) where {T <: CoordinateTransformations.AffineMap}
  LinearMap( convert(_Rot.QuatRotation, x) )
end

function convert(::Type{T}, x::SE3) where {T <: CoordinateTransformations.AffineMap}
  Translation(x.t...) ∘ convert(AffineMap{_Rot.QuatRotation{Float64}}, x.R)
end
function convert(::Type{SE3}, x::T) where {T <: CoordinateTransformations.AffineMap{_Rot.QuatRotation{Float64}}}
  SE3(x.translation[1:3], TransformUtils.Quaternion(x.linear.w, [x.linear.x,x.linear.y,x.linear.z]) )
end

function convert(::Type{SE3}, t::Tuple{Symbol, Vector{Float64}})
  if t[1]==:XYZqWXYZ
    return SE3(t[2][1:3],TransformUtils.Quaternion(t[2][4],t[2][5:7]))
  else
    error("Unknown conversion type $(t[1])")
  end
end

function convert(::Type{RangeAzimuthElevation}, val::Tuple{Symbol, Vector{Float64}})
  if val[1] == :rangeazimuth
    return RangeAzimuthElevation(val[2][1],val[2][2],nothing)
  elseif val[1] == :rangeazimuthelevation
    return RangeAzimuthElevation(val[2][1],val[2][2],val[2][3])
  else
    error("Unknown conversion from $(val[1]) to RangeAzimuthElevation.")
  end
end

# should be deprecated or indicated more clearly
lsrBR(a) = [a[2,:];a[1,:]]';


function \(s::SE3, wTr::CTs.Translation)
  bTr = s.R.R'*(wTr.v-s.t)
  Dtr = bTr
  range = norm(Dtr)
  azi = atan(Dtr[2], Dtr[1])
  elev = atan(Dtr[3], Dtr[1])
  RangeAzimuthElevation(range, azi, elev)
end



"""
    $SIGNATURES

Method to compare current and predicted estimate on a variable, developed for testing a new factor before adding to the factor graph.

Notes
- `fct` does not have to be in the factor graph -- likely used to test beforehand.
- function is useful for detecting if `multihypo` should be used.
- `approxConv` will project the full belief estimate through some factor but must already be in factor graph.

Example

```julia
# fg already exists containing :x7 and :l3
pp = Pose2Point2BearingRange(Normal(0,0.1),Normal(10,1.0))
# possible new measurement from :x7 to :l3
curr, pred = predictVariableByFactor(fg, :l3, pp, [:x7; :l3])
# example of naive user defined test on fit score
fitscore = minkld(curr, pred)
# `multihypo` can be used as option between existing or new variables
```

Related

approxConv
"""
function predictVariableByFactor( dfg::AbstractDFG,
                                  targetsym::Symbol,
                                  fct::Pose2Point2BearingRange,
                                  prevars::Vector{Symbol}  )
  #
  @assert targetsym in prevars
  curr = getBelief(dfg, targetsym)
  tfg = initfg()
  for var in prevars
    varnode = getVariable(dfg, var)
    addVariable!(tfg, var, getSofttype(varnode))
    if var != targetsym
      @assert isInitialized(varnode)
      initManual!(tfg,var,getBelief(varnode))
    end
  end
  addFactor!(tfg, prevars, fct, graphinit=false)
  fctsym = ls(tfg, targetsym)

  pts, infd = predictbelief(tfg, targetsym, fctsym)
  pred = manikde!(getManifold(getVariable(dfg, targetsym)), pts)
  # return current and predicted beliefs
  return curr, pred
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

function measureMeanDist(fg::AbstractDFG, a::AbstractString, b::AbstractString)
    #bearrang!(residual::Array{Float64,1}, Z::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
    res = zeros(2)
    A = getVal(fg,a)
    B = getVal(fg,b)
    Ax = Statistics.mean(vec(A[1,:]))
    Ay = Statistics.mean(vec(A[2,:]))
    Bx = Statistics.mean(vec(B[1,:]))
    By = Statistics.mean(vec(B[2,:]))
    dx = Bx - Ax
    dy = By - Ay
    b = atan(dy,dx)
    r = sqrt(dx^2 + dy^2)
    return r, b
end

function predictBodyBR(fg::AbstractDFG, a::Symbol, b::Symbol)
  res = zeros(2)
  A = getKDEMean(getBelief(getVariable(fg,a)))
  B = getKDEMean(getBelief(getVariable(fg,b)))
  Ax = A[1] # Statistics.mean(vec(A[1,:]))
  Ay = A[2] # Statistics.mean(vec(A[2,:]))
  Ath = getKDEMax(getBelief(getVariable(fg,a)))[3]
  Bx = B[1]
  By = B[2]
  wL = SE2([Bx;By;0.0])
  wBb = SE2([Ax;Ay;Ath])
  bL = se2vee((wBb \ Matrix{Float64}(LinearAlgebra.I, 3,3)) * wL)
  dx = bL[1] - 0.0
  dy = bL[2] - 0.0
  b = (atan(dy,dx))
  r = sqrt(dx^2 + dy^2)
  return b, r
end


"""
    $SIGNATURES

Return the last `number::Int` of poses according to `filterLabel::Regex`.

Notes
- Uses FIFO add history of variables to the distribued factor graph object as search index.
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

Create a new variable node and insert odometry constraint factor between
which will automatically increment latest pose symbol x<k+1> for new node new node and
constraint factor are returned as a tuple.
"""
function addOdoFG!(fg::G,
                   n::Symbol,
                   DX::Array{Float64,1},
                   cov::Array{Float64,2};
                   N::Int=0,
                   solvable::Int=1,
                   labels::Vector{<:AbstractString}=String[]  ) where G <: AbstractDFG
    #
    prev, X, nextn = getLastPose2D(fg)
    r,c = size(X)
    if N==0
      N = c
    end
    sig = diag(cov)
    XnextInit = zeros(r,c)
    # increases the number of particles based on the number of modes in the measurement Z
    for i in 1:c
        ent = [randn()*sig[1]; randn()*sig[2]; randn()*sig[3]]
        XnextInit[:,i] = addPose2Pose2(X[:,i], DX + ent)
    end

    v = addVariable!(fg, n, Pose2, N=N, solvable=solvable, tags=[labels;"POSE"])
    # v = addVariable!(fg, n, XnextInit, cov, N=N, solvable=solvable, tags=labels)
    pp = Pose2Pose2(MvNormal(DX, cov)) #[prev;v],
    f = addFactor!(fg, [prev;v], pp, solvable=solvable, graphinit=true )
    infor = inv(cov^2)
    # addOdoRemote(prev.index,v.index,DX,infor) # this is for remote factor graph ref parametric solution -- skipped internally by global flag variable
    return v, f
end

function addOdoFG!(fgl::G,
                   Z::Pose3Pose3;
                   N::Int=0,
                   solvable::Int=1,
                   labels::Vector{<:AbstractString}=String[]  ) where {G <: AbstractDFG}
  #
  vprev, X, nextn = getLastPoses(fgl)[1]
  vnext = addVariable!(fgl, nextn, Pose3, solvable=solvable, tags=labels)
  fact = addFactor!(fgl, [vprev;vnext], Z, graphinit=true)

  return vnext, fact

  # error("addOdoFG!( , ::Pose3Pose3, ) not currently usable, there were breaking changes. Work in Progress")
  # addOdoFG(fg, n, DX, cov, N=N, solvable=solvable, tags=labels)
end

"""
    $(SIGNATURES)

Create a new variable node and insert odometry constraint factor between
which will automatically increment latest pose symbol x<k+1> for new node new node and
constraint factor are returned as a tuple.

"""
function addOdoFG!(
        fgl::AbstractDFG,
        odo::Pose2Pose2;
        N::Int=0,
        solvable::Int=1,
        labels::Vector{<:AbstractString}=String[] ) # where {PP <: RoME.BetweenPoses}
    #
    vprev, X, nextn = getLastPose(fgl)
    if N==0
      N = size(X,2)
    end
    # vnext = addVariable!(fgl, nextn, X⊕odo, ones(1,1), N=N, solvable=solvable, tags=labels)
    vnext = addVariable!(fgl, nextn, Pose2, N=N, solvable=solvable, tags=labels)
    fact = addFactor!(fgl, [vprev;vnext], odo, graphinit=true)

    return vnext, fact
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




function malahanobisBR(measA, preA, cov::Array{Float64,2})
    # measure landmark with noise
    res = measA - preA
    mala2 = Union{}
    #Malahanobis distance
    if false
      lambda = cov \ Matrix{Float64}(LinearAlgebra.I, 2,2)
      mala2 = res' * lambda * res
    else
      mala2 = res' * (cov \ res)
    end

    mala = sqrt(mala2)
    return mala
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
      if length( DFG.getNeighbors(fg, id ) ) >= minnei
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
    if length( DFG.getNeighbors(fg, id) ) >= minnei
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

@deprecate getVariablesLabelsWithinRange(fg::AbstractDFG,regexKey::Regex=r"x";from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),minnei::Int=0) listVariablesLabelsWithinRange(fg,regexKey,from=from,to=to,minnei=minnei)

function get2DSampleMeans(fg::AbstractDFG,
                          regexKey::Regex=r"x";
                          from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),
                          minnei::Int=0)
  #
  X = Array{Float64,1}()
  Y = Array{Float64,1}()
  Th = Array{Float64,1}()
  LB = String[]

  vsyms = getVariablesLabelsWithinRange(fg, regexKey, from=from, to=to, minnei=minnei)

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
