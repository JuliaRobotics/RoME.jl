
## TODO: revamp these older functions to newer API in IIF/RoME/Caesar (Nov2018)
# standardize a common front-end design with feature tracking and real-time integration
# part of DFG discussion

function triggerPose(x, xprev, Tnow, Tprev,
                    distrule, timerule, yawrule)

  if norm(x[1:2]-xprev[1:2]) >= distrule
    @show Tnow, round(x,digits=2), round(xprev,digits=2)
    @show norm(x[1:2]-xprev[1:2])
    return 1
  end
  if abs(x[3]-xprev[3]) >= yawrule
    @show Tnow, round(x,digits=2), round(xprev,digits=2)
    @show abs(x[3]-xprev[3])
    return 2
  end
  if Tnow-Tprev > timerule
    return 3
  end
  return 0
end

mutable struct GenericInSituSystem{T}
  xprev::Array{Float64,1}
  x::Array{Float64,1}
  dOdo::Dict{Int,Array{Float64,1}}
  FeatAssc::Dict{Int, Dict{Int,Array{Float64,1}}}
  Tprev::Float64
  T0::Float64
  poseid::Int
  wTbk1::Array{Float64,2}
  bk1Tbk::Array{Float64,2}
  lstlaseridx::Int
  trackers::Dict{Int,T}
end

const InSituSystem = GenericInSituSystem{Feature}

function makeInSituSys(x::Array{Float64,1}, bfts0::Array{Float64,2})
  dOdo = Dict{Int,Array{Float64,1}}()
  FeatAssc = Dict{Int, Dict{Int,Array{Float64,1}}}()
  Tprev = 0.0
  T0 = 0.0
  wTbk1 = SE2(x)
  bk1Tbk = SE2(zeros(3))
  poseid = 1
  dOdo[poseid] = [x[1];x[2];x[3];T0;0]
  lstlaseridx = 1
  trackers = initTrackersFrom(bfts0)
  return InSituSystem(
  x,
  x,
  dOdo,
  FeatAssc,
  Tprev,
  T0,
  poseid,
  wTbk1,
  bk1Tbk,
  lstlaseridx,
  trackers
  )
end


function makeGenericInSituSys(x::Array{Float64,1})
  dOdo = Dict{Int,Array{Float64,1}}()
  FeatAssc = Dict{Int, Dict{Int,Array{Float64,1}}}()
  Tprev = 0.0
  T0 = 0.0
  wTbk1 = SE2(x)
  bk1Tbk = SE2(zeros(3))
  poseid = 1
  dOdo[poseid] = [x[1];x[2];x[3];T0;0]
  lstlaseridx = 1
  trackers = Dict{Int,Any}()
  return GenericInSituSystem{Any}(
  x,
  x,
  dOdo,
  FeatAssc,
  Tprev,
  T0,
  poseid,
  wTbk1,
  bk1Tbk,
  lstlaseridx,
  trackers
  )
end

# doesn't work in full general case yet, still requires reset to zero at each new pose
function poseTrigAndAdd!(instSys::GenericInSituSystem, Ts::Float64,
                        distrule::Float64, timerule::Float64, yawrule::Float64;
                        xprev=zeros(3), auxtrig::Bool=false)
  rule = triggerPose(instSys.x, xprev, Ts, instSys.Tprev, distrule, timerule, yawrule)
  if rule != 0 || auxtrig
    instSys.bk1Tbk = SE2(instSys.x)
    n = [instSys.x;[Ts;rule]]
    instSys.poseid +=1
    instSys.dOdo[instSys.poseid] = n
    instSys.wTbk1 = instSys.wTbk1*instSys.bk1Tbk
    instSys.Tprev = Ts
    instSys.x[1] = 0.0; instSys.x[2]=0.0; instSys.x[3]=0.0
    return true
  end
  return false
end

function poseTrigAndAdd!(instSys::InSituSystem, Ts::Float64,
                        distrule::Float64, timerule::Float64, yawrule::Float64;
                        xprev=zeros(3), auxtrig::Bool=false)

  error(" Yeah dehann, do it properly")


end


function processTreeTrackersUpdates!(instSys::InSituSystem, lsrFeats::Dict{Int,LaserFeatures},
                                    Ts::Float64, b1Dxb::Array{Float64,1}, DBG=Union{}; dbgflag=true)
  propAllTrackers!(instSys.trackers, b1Dxb, [0.05;0.05;0.004])
  newlsridx, Ta = getFeatsAtT(lsrFeats, Ts, prev=instSys.lstlaseridx)
  if newlsridx != instSys.lstlaseridx
    # we have a new laser scan available and should update the trackers
    instSys.lstlaseridx = newlsridx
    bfts = lsrFeats[instSys.lstlaseridx].feats
    hardassc = assocMeasWFeats!(instSys.trackers, bfts)
    measUpdateTrackers!(instSys.trackers, hardassc, [0.5;0.05]) # TODO, noise here is not in polar coords!! something wrong
    if dbgflag
      DBG[instSys.lstlaseridx] = deepcopy(instSys.trackers)
    end
  end

  nothing
end

function advOdoByRules(DRS::Array{Float64,2}, lsrFeats::Dict{Int,LaserFeatures};
                        distrule=20.0, timerule=30.0, yawrule=pi/3.0, trkfeats=true)
  # DBG = Dict{Int,Dict{Int,Array{BallTreeDensity,1}}}()
  # DBG = Dict{Int,Dict{Int,Feature}}()
  DBG = Union{}

  bfts0 = lsrFeats[1].feats
  instSys = makeInSituSys(zeros(3), bfts0)
  fdict = Dict{Int,Array{Float64,1}}()
  for f in instSys.trackers
    # if length(f[2].lastz) < 3
    #   error("HEY!! $(f[1])")
    # end
    fdict[f[2].id] = f[2].lastz
  end
  instSys.FeatAssc[instSys.poseid] = fdict
  for i in 1:size(DRS,1)
    # do odo
    dt = DRS[i,1]-instSys.T0
    whlspd, strang = compensateRawDRS(vec(DRS[i,:]))
    bTbm = SE2(instSys.x)
    instSys.x=uteOdomEasy(instSys.x, whlspd, strang, dt)
    bTbp = SE2(instSys.x)

    # propagate, and maybe measure update, feature trackers
    if trkfeats
      #@show keys(instSys.trackers)
      bmTbp = vec(se2vee(inv(bTbm)*bTbp))
      @time processTreeTrackersUpdates!(instSys, lsrFeats, DRS[i,1], bmTbp, DBG,dbgflag=false)
    end
    # subsample motion tracks for constructing the factor graph SLAM solution
    if poseTrigAndAdd!(instSys, DRS[i,1], distrule, timerule, yawrule)
      fdict = Dict{Int,Array{Float64,1}}()
      for f in instSys.trackers
        mpt = vec(Base.mean(getPoints(f[2].bel),2))
        r,b = c2p(mpt)
        fdict[f[2].id] = [r;b;f[2].lastz[3]]
      end
      instSys.FeatAssc[instSys.poseid] = fdict
    end
    # advance previous time point
    instSys.T0 = DRS[i,1]
  end
  return instSys.dOdo, instSys.FeatAssc, DBG
end




# function loadVicPrkDataset(filename::AbstractString="datasets/VicPrk.jld")
#   DRS,GPS,LsrFeats,d,f = jldopen(filename, "r") do file
#     read(file, "DRS")
#     read(file, "GPS")
#     read(file, "LsrFeats")
#     read(file, "d")
#     read(file, "f")
#   end
#   return DRS,GPS,LsrFeats,d,f
# end
