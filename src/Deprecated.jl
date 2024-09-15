
##==============================================================================
## Legacy, remove once AMP #41 is resolved
##==============================================================================


Base.convert(
  ::Type{<:Tuple}, 
  ::IIF.InstanceType{typeof(getManifold(RotVelPos))}
) = (
  :Circular,:Circular,:Circular,
  :Euclid,:Euclid,:Euclid,
  :Euclid,:Euclid,:Euclid
)


##==============================================================================
## Legacy, remove some time after RoME v0.22
##==============================================================================

# @deprecate homographyToCoordinates(w...;kw...) homography_to_coordinates(w...;kw...)

# export measureMeanDist
# function measureMeanDist(fg::AbstractDFG, a::AbstractString, b::AbstractString)
#   @error "RoME.measureMeanDist is obsolete"
#   #bearrang!(residual::Array{Float64,1}, Z::Array{Float64,1}, X::Array{Float64,1}, L::Array{Float64,1})
#   res = zeros(2)
#   A = getVal(fg,a)
#   B = getVal(fg,b)
#   Ax = Statistics.mean(vec(A[1,:]))
#   Ay = Statistics.mean(vec(A[2,:]))
#   Bx = Statistics.mean(vec(B[1,:]))
#   By = Statistics.mean(vec(B[2,:]))
#   dx = Bx - Ax
#   dy = By - Ay
#   b = atan(dy,dx)
#   r = sqrt(dx^2 + dy^2)
#   return r, b
# end

# # should be deprecated or indicated more clearly
# @deprecate lsrBR(a) [a[2,:];a[1,:]]'

# # import AMP: _makeVectorManifold
# # AMP._makeVectorManifold(::M, prr::ProductRepr) where {M <: typeof(BearingRange_Manifold)} = coords(M, prr)



##==============================================================================
## Remove as part of Manifolds.jl consolidation, #244
##==============================================================================

# export veePose3, veePose

# function veePose3(s::SE3)
#   TransformUtils.veeEuler(s)
# end
# function veePose(s::SE3)
#   TransformUtils.veeEuler(s)
# end


# legacy support, will be deprecated
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{typeof(AMP.SE2_Manifold)}) = (:Euclid, :Euclid, :Circular)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{typeof(SE2E2_Manifold)}) = (:Euclid,:Euclid,:Circular,:Euclid,:Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{typeof(BearingRange_Manifold)}) = (:Circular,:Euclid)
Base.convert(::Type{<:Tuple}, ::IIF.InstanceType{typeof(getManifold(DynPose2))}) = (:Euclid,:Euclid,:Circular,:Euclid,:Euclid)


Base.convert(
  ::Type{<:Tuple},
  ::IIF.InstanceType{typeof(Manifolds.ProductGroup(ProductManifold(SpecialEuclidean(2), TranslationGroup(2)), LeftInvariantRepresentation()))}
) = (:Euclid,:Euclid,:Circular,:Euclid,:Euclid)


# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Point2Point2}) = AMP.Euclid2
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Pose2Point2}) = AMP.Euclid2
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Pose2Point2Bearing}) = AMP.Euclid
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Point2Point2Range}) = AMP.Euclid
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Pose2Point2Range}) = AMP.Euclid
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Pose2Point2BearingRange}) = AMP.Euclid2
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Pose2Pose2}) = AMP.SE2_Manifold
# Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{Pose3Pose3}) = AMP.SE3_Manifold
Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{DynPoint2DynPoint2}) = AMP.Euclid4
Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{DynPose2DynPose2}) = begin @warn("FIXME"); SE2E2_Manifold end
Base.convert(::Type{<:ManifoldsBase.AbstractManifold}, ::IIF.InstanceType{VelPose2VelPose2}) = begin @warn("FIXME"); SE2E2_Manifold end



##==============================================================================
## OLD Victory Park Example code, don't delete until replacements are coded
##==============================================================================



# function projNewLandmPoints(vps::Graphs.ExVertex, br::Array{Float64,1}, cov::Array{Float64,2})
#   # TODO -- convert to use Distributions and common projection function
#   Xps = getVal(vps)
#   lmPts = zeros(2,size(Xps,2))
#   for i in 1:size(Xps,2)
#       ent = [cov[1,1]*randn(); cov[2,2]*randn()]
#       init = vec(Xps[1:2,i])+randn(2)
#       lmPts[:,i] = solveLandm(br + ent, vec(Xps[:,i]), init)
#   end
#   return lmPts
# end



# function calcIntersectVols( fgl::G, predLm::BallTreeDensity;
#                             currage=0, maxdeltaage=Inf) where G <: AbstractDFG
#     # TODO upgrade to using MMD test
#     # all landmarks of interest
#     ll = ls(fgl, r"l\d")
#     # output result
#     rr = Dict{String, RemoteRef}()
#     fetchlist = String[]
#     iv = Dict{String, Float64}()
#     for l in ll
#       pvlm = getVariable(fgl,l)
#       # TODO -- can be improved via query in DB case
#       if currage - pvlm.attributes["maxage"] < maxdeltaage
#         p = getBelief(fgl, l)
#         rr[l] = remotecall(uppA(), intersIntgAppxIS, p,predLm)
#         push!(fetchlist, l)
#       else
#         println("calcIntersectVols -- ignoring $(l) because maxdeltaage exceeded")
#         iv[l] = 0
#       end
#     end
#     max = 0
#     maxl = String("")
#     for l in fetchlist #ll
#       # p = getBelief(fgl, l)
#       # tv = intersIntgAppxIS(p,predLm)
#       # iv[l] = tv
#       tv = fetch(rr[l])
#       iv[l] = tv
#       if max < tv  max = tv; maxl = l; end
#     end
#     return iv, maxl
# end

# function maxIvWithoutID(ivs::Dict{String, Float64}, l::T) where {T <: AbstractString}
#   max = 0
#   maxl = String("")
#   for i in ivs
#     if max < i[2] && i[1] != l;  max = i[2]; maxl = i[1]; end
#   end
#   return maxl
# end

# # binary tests to distinguish how to automatically add a landmark to the existing factor graph
# function doAutoEvalTests(fgl::G, ivs::Dict{T, Float64}, maxl::T, lmid::Int, lmindx::Int) where {G <: AbstractDFG, T <: AbstractString}
#   maxAnyval = maxl != String("") ? ivs[maxl] : 0.0
#   # maxid = fgl.IDs[maxl]
#   lmidSugg = lmid != -1 # a landmark ID has been suggested
#   maxAnyExists = maxAnyval > 0.03 # there is notable intersection with previous landm
#   lmIDExists = haskey(fgl.v, lmid) # suggested lmid already in fgl
#   newlmindx = lmindx
#   if lmIDExists
#     lmSuggLbl = String(getVert(fgl,lmid).label) # TODO -- wasteful
#   else
#     newlmindx = lmindx + 1
#     lmSuggLbl = String(string('l',newlmindx))
#   end
#   maxl2 = lmIDExists ? maxIvWithoutID(ivs, lmSuggLbl) : String("")
#   maxl2Exists = lmIDExists ? (maxl2 != "" ? ivs[maxl2] > 0.03 : false) : false # there is notable intersection with previous landm
#   intgLmIDExists = lmIDExists ? ivs[lmSuggLbl] > 0.03 : false

#   return lmidSugg, maxAnyExists, maxl2Exists, maxl2, lmIDExists, intgLmIDExists, lmSuggLbl, newlmindx
# end


# function evalAutoCases!(fgl::G, lmid::Int, ivs::Dict{T, Float64}, maxl::T,
#                         pose::T, lmPts::Array{Float64,2}, br::Array{Float64,1}, cov::Array{Float64,2}, lmindx::Int;
#                         N::Int=100, solvable::Int=1 ) where {G <: AbstractDFG, T <: AbstractString}
#   lmidSugg, maxAnyExists, maxl2Exists, maxl2, lmIDExists, intgLmIDExists, lmSuggLbl, newlmindx = doAutoEvalTests(fgl,ivs,maxl,lmid, lmindx)

#   println("evalAutoCases -- found=$(lmidSugg), $(maxAnyExists), $(maxl2Exists), $(lmIDExists), $(intgLmIDExists)")

#   vlm = Union{}; fbr = Union{};
#   if (!lmidSugg && !maxAnyExists)
#     #new landmark and UniBR constraint
#     v,L,lm = getLastLandm2D(fgl)
#     vlm = newLandm!(fgl, lm, lmPts, cov, N=N,solvable=solvable)
#     fbr = addBRFG!(fgl, pose, lm, br, cov, solvable=solvable)
#   elseif !lmidSugg && maxAnyExists
#     # add UniBR to best match maxl
#     vlm = getVariable(fgl,maxl)
#     fbr = addBRFG!(fgl, pose, maxl, br, cov, solvable=solvable)
#   elseif lmidSugg && !maxl2Exists && !lmIDExists
#     #add new landmark and add UniBR to suggested lmid
#     vlm = newLandm!(fgl, lmSuggLbl, lmPts, cov, N=N, solvable=solvable)
#     fbr = addBRFG!(fgl, pose, lmSuggLbl, br, cov, solvable=solvable)
#   elseif lmidSugg && !maxl2Exists && lmIDExists && intgLmIDExists
#     # doesn't self intesect with existing lmid, add UniBR to lmid
#     vlm = getVariable(fgl, lmid)
#     fbr = addBRFG!(fgl, pose, lmSuggLbl, br, cov, solvable=solvable)
#   elseif lmidSugg && maxl2Exists && !lmIDExists
#     # add new landmark and add MMBR to both maxl and lmid
#     vlm = newLandm!(fgl, lmSuggLbl, lmPts, cov, N=N, solvable=solvable)
#     addMMBRFG!(fgl, pose, [maxl2;lmSuggLbl], br, cov, solvable=solvable)
#   elseif lmidSugg && maxl2Exists && lmIDExists && intgLmIDExists
#     # obvious case, add MMBR to both maxl and lmid. Double intersect might be the same thing
#     println("evalAutoCases! -- obvious case is happening")
#     addMMBRFG!(fgl, pose, [maxl2;lmSuggLbl], br, cov, solvable=solvable)
#     vlm = getVariable(fgl,lmSuggLbl)
#   elseif lmidSugg && maxl2Exists && lmIDExists && !intgLmIDExists
#     # odd case, does not intersect with suggestion, but does with some previous landm
#     # add MMBR
#     @warn "evalAutoCases! -- no self intersect with suggested $(lmSuggLbl) detected"
#     addMMBRFG!(fgl, pose, [maxl;lmSuggLbl], br, cov, solvable=solvable)
#     vlm = getVariable(fgl,lmSuggLbl)
#   elseif lmidSugg && !maxl2Exists && lmIDExists && !intgLmIDExists
#   #   # landm exists but no intersection with existing or suggested lmid
#   #   # may suggest some error
#     @warn "evalAutoCases! -- no intersect with suggested $(lmSuggLbl) or map detected, adding  new landmark MM constraint incase"
#     v,L,lm = getLastLandm2D(fgl)
#     vlm = newLandm!(fgl, lm, lmPts, cov, N=N, solvable=solvable)
#     addMMBRFG!(fgl, pose, [lm; lmSuggLbl], br, cov, solvable=solvable)
#   else
#     error("evalAutoCases! -- unknown case encountered, can reduce to this error to a warning and ignore user request")
#   end

#   return vlm, fbr, newlmindx
# end

# function addAutoLandmBR!(fgl::G,
#                          pose::T,
#                          lmid::Int,
#                          br::Array{Float64,1},
#                          cov::Array{Float64,2},
#                          lmindx::Int;
#                          N::Int=100,
#                          solvable::Int=1  ) where {G <: AbstractDFG, T <: AbstractString}
#   #
#   vps = getVariable(fgl, pose)
#   lmPts = projNewLandmPoints(vps, br, cov)
#   lmkde = kde!(lmPts)
#   currage = parse(Int, pose[2:end])
#   ivs, maxl = calcIntersectVols(fgl, lmkde, currage=currage,maxdeltaage=10)

#   # There are 8 cases of interest
#   vlm, fbr, newlmindx = evalAutoCases!(fgl, lmid, ivs, maxl,pose,lmPts, br,cov,lmindx,N=N,solvable=solvable)

#   return vlm, fbr, newlmindx
# end


# function newLandm!(fg::AbstractDFG, lm::T, wPos::Array{Float64,2}, sig::Array{Float64,2};
#                   N::Int=100, solvable::Int=1, labels::Vector{T}=String[]) where {T <: AbstractString}

#     vert=addVariable!(fg, Symbol(lm), Point2, N=N, solvable=solvable, tags=union(["LANDMARK";], labels))
#     # TODO -- need to confirm this function is updating the correct memory location. v should be pointing into graph
#     # vert=addVariable!(fg, Symbol(lm), wPos, sig, N=N, solvable=solvable, tags=labels)

#     vert.attributes["age"] = 0
#     vert.attributes["maxage"] = 0
#     vert.attributes["numposes"] = 0
#     updateFullVert!(fg, vert)

#     println("newLandm! -- added $(lm)")
#     return vert
# end

# function addBRFG!(fg::G,
#                   pose::T,
#                   lm::T,
#                   br::Array{Float64,1},
#                   cov::Array{Float64,2};
#                   solvable::Int=1  ) where {G <: AbstractDFG, T <: AbstractString}
#   #
#   vps = getVert(fg,pose)
#   vlm = getVert(fg,lm)
#   testlbl = vps.label*vlm.label
#   for nei in getOutNeighbors(fg, vlm)
#     if nei.label == testlbl
#       # TODO -- makes function call brittle
#       @warn "We already have $(testlbl), skipping this constraint"
#       return nothing
#     end
#   end
#   @show keys(vlm.attributes)
#   np = vlm.attributes["numposes"]
#   la = vlm.attributes["age"]
#   nage = parse(Int,pose[2:end])
#   vlm.attributes["numposes"] = np+1
#   vlm.attributes["age"] = ((la*np)+nage)/(np+1)
#   vlm.attributes["maxage"] = nage
#   updateFullVert!(fg, vlm)


#   pbr = Pose2Point2BearingRange(Normal(br[1], cov[1,1]), Normal(br[2],  cov[2,2]))  #{Normal, Normal}
#   @show vps, vlm
#   f = addFactor!(fg, [vps;vlm], pbr, solvable=solvable, graphinit=true ) #[vps;vlm],

#   # only used for max likelihood unimodal tests.
#   u, P = pol2cart(br[[2;1]], diag(cov))
#   infor = inv(P^2)
#   # addLandmMeasRemote(vps.index,vlm.index,u,infor) # for iSAM1 remote solution as reference
#   return f
# end

# function addMMBRFG!(fg::G,
#                     syms::Array{Symbol,1}, br::Array{Float64,1},
#                     cov::Array{Float64,2}; w::Vector{Float64}=Float64[0.5;0.5],
#                     solvable::Int=1) where G <: AbstractDFG
#     #
#     # vps = getVert(fg,pose)
#     # vlm1 = getVert(fg,lm[1])
#     # vlm2 = getVert(fg,lm[2])

#     pbr = Pose2Point2BearingRange(Normal(br[1],cov[1,1]),  Normal(br[2],cov[2,2]))
#     syms = Symbol.([pose;lm...])
#     f = addFactor!(fg, syms, pbr, multihypo=[1.0; w...], solvable=solvable, graphinit=true )
#     return f
# end

# function projNewLandm!(fg::G,
#                        pose::T,
#                        lm::T,
#                        br::Array{Float64,1},
#                        cov::Array{Float64,2};
#                        addfactor=true,
#                        N::Int=100,
#                        solvable::Int=1,
#                        labels::Vector{T}=String[]  ) where {G <: AbstractDFG, T <: AbstractString}
#     #
#     vps = getVariable(fg, pose)

#     lmPts = projNewLandmPoints(vps, br, cov)
#     vlm = newLandm!(fg, lm, lmPts, cov, N=N, solvable=solvable, tags=labels) # cov should not be required here
#     if addfactor
#       fbr = addBRFG!(fg, pose, lm, br, cov, solvable=solvable)
#       return vlm, fbr
#     end
#     return vlm
# end




##==============================================================================
## OLD Victory Park Example code, don't delete until replacements are coded
##==============================================================================



# # optional tools
# using Requires

# function __init__()
#   # combining neural networks natively into the non-Gaussian  factor graph object
#   # @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin
#   #   # include("factors/flux/models/Pose2OdoNN_01.jl") # until a better way is found to deserialize
#   #   # include("factors/flux/MixtureFluxPose2Pose2.jl")
#   # end

#   # Scalar field specifics
  
#   # @require ImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534" begin
#   #   @require ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19" include("services/RequiresImages.jl")
#   # end
#   # Images="916415d5-f1e6-5110-898d-aaa5f9f070e0" 

#   # @require Interpolations="a98d9a8b-a2ab-59e6-89dd-64a1c18fca59" begin 
#   #   include("services/ScalarFieldsInterpolations.jl")
#   # end
# end