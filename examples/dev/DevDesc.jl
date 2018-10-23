include("../hopfield.jl")
using Gadfly, Colors

function circLayer(x::Vector{Float64};r::Float64=1.0,c="red")
  line = 0:0.05:2pi
  N = length(line)
  X, Y = zeros(N), zeros(N)
  for i in 1:N
    X[i] = x[1] + r*cos(line[i])
    Y[i] = x[2] + r*sin(line[i])
  end
  return Gadfly.layer(x=X,y=Y, Geom.path(), Gadfly.Theme(default_color=parse(Colorant,c)))
end

function drawLandmPtsLayer(pts::Dict{String, Vector{Float64}}; c::String="blue")
  N = length(pts)
  xy = zeros(2,N)
  lbls = String[]
  i = 0
  lbl = String[]
  for p in pts
    i += 1
    push!(lbls, p[1])
    [xy[j,i] = p[2][j] for j in 1:2]
    push!(lbl, p[1])
  end
  return layer(x=xy[1,:],y=xy[2,:], label=lbl, Geom.point, Geom.label, Gadfly.Theme(default_color=parse(Colorant,c)))
end

incirc(dx,dy,r) = dx^2 + dy^2 <= r^2

function findInRadius(pts::Dict{String, Vector{Float64}}, from::String; r::Float64=1.0)
  rd = Dict{String, Vector{Float64}}()
  xy = pts[from]
  rd[from] = pts[from]
  for p in pts
    if incirc(xy[1]-p[2][1], xy[2]-p[2][2], r)
      rd[p[1]] = p[2]
    end
  end
  return rd
end

function drawLineLayer(from::Vector{Float64}, to::Vector{Float64}; c="magenta")
  return Gadfly.layer(x=[from[1];to[1]],y=[from[2];to[2]], Geom.path(), Gadfly.Theme(default_color=parse(Colorant,c)))
end

function drawDictLinesLayers(pl, pts::Dict{String, Vector{Float64}}, from::String)
  x0 = pts[from][1:2]
  for p in pts
    p[1] == from ? continue : nothing
    l = drawLineLayer(x0, p[2][1:2])
    push!(pl.layers, l[1])
  end
  return pl
end

calcang(x1::Vector{Float64},x2::Vector{Float64}) = atan(x2[2]-x1[2], x2[1]-x1[1])
calcrange(x1::Vector{Float64},x2::Vector{Float64}) = norm(x1-x2)

function orderedBinIndex(x::Float64; start::Float64=0.0, stop::Float64=1.0,N::Int=3)
  refSpace= range(start, stop=stop, length=N+1)
  if N == 1
    return 1
  end
  for i in 2:(N+1)
    if refSpace[i-1] <= x < refSpace[i]
      return i-1
    end
  end
  return N
end

function constructShapeContexctMat(range::Dict{String, Float64}, ang::Dict{String, Float64}, me::String;
                                   angBins::Int=4, rangeBins::Int=3, r::Float64=1.0)
  sc = zeros(rangeBins, angBins)
  for a in ang
    if String(a[1]) == me  continue; end
    # @show a[1], a[2], range[a[1]]
    aIdx = orderedBinIndex(a[2], start=Float64(-pi), stop=Float64(pi), N=angBins)
    rIdx = orderedBinIndex(range[a[1]], start=0.0, stop=r, N=rangeBins)
    # @show rIdx, aIdx
    sc[rIdx, aIdx] = sc[rIdx, aIdx] + 1
  end
  return sc
end

function myShapeContext(pts::Dict{String, Vector{Float64}}, me::String, refpt::String;
                        r::Float64=1.0, angBins::Int=4, rangeBins::Int=3)
  rd = findInRadius(pts, me, r=r)
  ang = Dict{String, Float64}()
  rang = Dict{String, Float64}()
  # scVec = Dict{String}
  refang = calcang(rd[me], pts[refpt])
  for p in rd
    ang[p[1]] = calcang(rd[me], p[2]) - refang
    rang[p[1]] = calcrange(rd[me], p[2])
  end
  # @show ang
  # @show rang
  mat = constructShapeContexctMat(rang, ang, me, r=r, angBins=angBins, rangeBins=rangeBins)
  return mat
end


function calcShapeContext(pts::Dict{String, Vector{Float64}}, refpt::String;
                          r::Float64=1.0, angBins::Int=4, rangeBins::Int=3)
  rd = findInRadius(pts,refpt,r=r)
  mat = zeros(rangeBins, angBins)
  for doone in keys(rd)
    if String(doone) == refpt  continue; end
    # @show doone
    mat = mat + myShapeContext(pts, doone, refpt, r=r, angBins=angBins, rangeBins=rangeBins)
  end
  return mat
end

function shapeContMatToVec(mat::Array{Float64,2}; duplicity::Int=2)
  v = [map(Int,mat[:].==1)'; map(Int,mat[:].>=2)']#; map(Int,mat[:].>2)']
  # v = map(Int,mat[:].>=1)'
  r,c = size(v)
  V = reshape(v,floor(Int,r*c),1)
  return sign(V[:]-0.5)
end

function shapeContPattern(pts::Dict{String, Vector{Float64}}, refpt::String;
                          r::Float64=1.0, angBins::Int=4, rangeBins::Int=3)
    mat = calcShapeContext(pts, refpt, r=r, angBins=angBins, rangeBins=rangeBins)
    return shapeContMatToVec(mat)
end




# we want a way to select the top 3 descriptors from current list of landmarks
# parameter
function bestDescriptors(d::Dict{String, Vector{Float64}};
                          r::Float64=1.0,Mf::Int=4,angBins::Int=4,rangeBins::Int=3,k::Int=3)

  allDescr = Dict{String, Vector{Int}}()
  skey = collect(keys(d))
  for lm in skey
    if length(findInRadius(d,lm,r=r)) > Mf
      allDescr[lm] = shapeContPattern(d, lm, r=r, angBins=angBins, rangeBins=rangeBins)
    end
  end
  @show ks = collect(keys(allDescr))

  #construct matrix of hamming distances
  n = length(ks)
  H = 99999*Matrix{Int}(LinearAlgebra.I, n, n)
  for i in 1:(n-1), j = (i+1):n
    H[i,j] = hamming(allDescr[ks[i]],allDescr[ks[j]])
  end
  H = H + H'

  # sort for max of min distances
  Hmin = vec(minimum(H,2))
  p = sortperm(Hmin,rev=true)
  # @show Hmin[p]
  # @show ks[p]

  # take top half and redo
  H1a = H[:,p]
  H1b = H1a[p,:]

  k = n >= k ? k : n
  # @show nn = k #floor(Int, n/2+0.5)
  H1 = H1b[1:k, 1:k]
  @show kss = (ks[p])[1:k]

  H1min = vec(minimum(H1,2))
  p1 = sortperm(H1min,rev=true)
  # @show H1min[p1]
  # @show kss[p1]
  bestlbls = kss[p1]#[1:k]]
  bestDescr = Dict{String, Vector{Int}}()
  for lb in bestlbls
    bestDescr[lb] = allDescr[lb]
  end
  return bestlbls, bestDescr
end


function findMatches(dd::Dict{String, Vector{Float64}}, Patt::Array{Float64,2}, bLbls::Vector{String};
  R::Float64=1.0,angBins::Int=4,rangeBins::Int=3, Mf::Int=4)
  SV = Dict{String, Vector{Int}}()
  Psv = Dict{String, Int}()
  for i in 1:length(dd)
    l = length(findInRadius(dd,"l$(i)",r=R))
    if l > Mf
      val = shapeContPattern(dd, "l$(i)", r=R, angBins=angBins, rangeBins=rangeBins)
      SV[String("l$(i)")] = val
      Psv[String("l$(i)")] = searchHopfield(Patt, val, maxiter=4)
    end
  end

  # display matching results
  mtchs = Dict{Int, Any}()
  fm = Dict{Int, Any}()
  # @show bLbls
  i = 0
  ord = Int[]
  for psv in Psv
    if psv[2] != -1
      i+=1
      hd = hamming(SV[psv[1]], Patt[:,psv[2]])
      mtchs[i] = (psv[1], bLbls[psv[2]], hd)
      push!(ord, hd)
    end
  end
  p = sortperm(ord)
  i = 0
  for j in p
    i+=1
    fm[i] = mtchs[j]
  end
  return fm
end


function convertToMMLst(fm::Dict{Int,Any}, distLim::Int)
  mmlst = Dict{String, Array{Tuple{String, Int},1}}()
  for i in 1:length(fm)
    if fm[i][3] <= distLim
      if haskey(mmlst,fm[i][2])
          push!(mmlst[fm[i][2]], (fm[i][1],fm[i][3]) )
      else
          mmlst[fm[i][2]] =  [ (fm[i][1],fm[i][3])  ]
      end
    end
  end
  return mmlst
end


# start some practicatical tests
if false
d = Dict{String, Vector{Float64}}()
for i in 1:9
  d[String("l$(i)")] = rand(2)
end
for i in 10:20
  d[String("l$(i)")] = (3.0.*rand(2))-1
end

# draw points
pl=plot(drawLandmPtsLayer(d))

R = 0.4
Mf = 4
angBins=6
bLbls, bestDescr = bestDescriptors(d, r=R, Mf=Mf, k=4 )

# draw so we can see
for lm in bLbls
  push!(pl,circLayer(d[lm],r=R)[1])
end
pl

# Training pattern vectors
# @show Plbl = [bLbls[1];bLbls[2];bLbls[3]]#; bLbls[4]]
@show bLbls
v1 = bestDescr[bLbls[1]]
v2 = bestDescr[bLbls[2]]
v3 = bestDescr[bLbls[3]]
v4 = bestDescr[bLbls[4]]

P = [v1';v2';v3';v4']'
P = map(Float64,P)

hamming(v1,v2)
hamming(v1,v3)
hamming(v2,v3)

hamming(v1,v4)
hamming(v2,v4)
hamming(v3,v4)



# how well do these patterns work

dd = deepcopy(d)
@show length(d)
for i in 1:length(d)
  dd[String("l$(i)")] = dd[String("l$(i)")] + 0.0*randn(2)
end
for i in length(d):floor(Int,length(d)*1.0)
  dd[String("l$(i)")] = (3.0.*rand(2))-1
end

findMatches(dd, P, bLbls, R=R, angBins=angBins, Mf=Mf)

l2 = drawLandmPtsLayer(dd,c="green")
push!(pl.layers, l2[1])
pl

# plot(drawLandmPtsLayer(dd))



# random drawing
refpt = "l8"
rd = findInRadius(d,refpt,r=R)
pl=plot(drawLandmPtsLayer(d),
     circLayer(d[refpt],r=R))
drawDictLinesLayers(pl, rd, refpt)


# sanity check

val1 = shapeContPattern(d, "l3", r=R, angBins=angBins)
val2 = shapeContPattern(dd, "l3", r=R, angBins=angBins)
@show hamming(val1,val2)

end

# some training and recognition tests
# ks = collect(keys(rd))
# doone = ks[1] != "l1" ? ks[1] : ks[2]
# doone="l6"
# push!(pl,circLayer(rd[doone],r=0.4)[1])
#
# mat = zeros(3, 4)
# for doone in ["l6";"l8";"l10"]
#   @show doone
#   @show mat = mat + myShapeContext(d, doone, "l1", r=0.4)
# end
#
# v1 = shapeContPattern(d, "l2", r=0.4)
# v2 = shapeContPattern(d, "l4", r=0.4)
# v3 = shapeContPattern(d, "l7", r=0.4)
# v4 = shapeContPattern(d, "l10", r=0.4)
#
# P = [v1';v2';v3']'
#
# a1 = deepcopy(v10)
# a1[1] = -1*v10[1]
# a1[10] = -1*v10[10]
# @show searchHopfield(P, a1)




# trying sorting

# @show A = rand(4,4)
#
# A = 16*Matrix{Float64}(LinearAlgebra.I, 4,4)
# for i in 1:3, j in (i+1):4
# A[i,j] = i*j
# end
# A = A + A'
# @show A
#
# [A[i,i]=24.0 for i in 1:4]
# @show Amin = vec(minimum(A,2))
# @show p = sortperm(Amin,rev=true)
# A1 = A[:,p]
# @show A2 = A1[p,:]
