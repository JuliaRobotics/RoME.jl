global DISABLESTBPRTLINES = false

function togglePrtStbLines()
  global DISABLESTBPRTLINES
  DISABLESTBPRTLINES = !DISABLESTBPRTLINES
end

function plotLsrScanFeats(br::Array{Float64,2})
  Cart = zeros(size(br))
  Cart[:,1] = br[:,2].*cos(br[:,1])
  Cart[:,2] = br[:,2].*sin(br[:,1])
  plot(x=Cart[:,1],y=Cart[:,2],Geom.point,
  Guide.xticks(ticks=collect(-60:10:60)),
  Guide.yticks(ticks=collect(0:10:80)))
end

function drawFeatTrackers(trkrs::Dict{Int,Feature}, bfts::Array{Float64,2})
  musX = Float64[]
  varX = Float64[]
  musY = Float64[]
  varY = Float64[]
  allPtsX = Float64[]
  allPtsY = Float64[]

  for ftr in trkrs
    pts = getPoints(ftr[2].bel)
    allPtsX = [allPtsX; vec(pts[1,:])]
    allPtsY = [allPtsY; vec(pts[2,:])]

    push!(musX, Base.mean(vec(pts[1,:])))
    push!(varX, Base.std(vec(pts[1,:])))
    push!(musY, Base.mean(vec(pts[2,:])))
    push!(varY, Base.std(vec(pts[2,:])))
  end

  X = Float64[]
  Y = Float64[]

  if size(bfts,2) > 0
    if bfts[1,1] != 0.0 && bfts[2,1] != 0.0 && bfts[3,1] != 0.0
      for i in 1:size(bfts,2)
          u, R = p2c(vec(bfts[:,i]))
          push!(X, u[1])
          push!(Y, u[2])
      end
    end
  end

  # Guide.yticks(ticks=collect(-60:10:60)),
  # Guide.xticks(ticks=collect(0:10:80))
  p = plot(layer(x=musX, y=musY, Geom.point, Theme(default_color=colorant"red")),
  layer(x=allPtsX, y=allPtsY, Geom.histogram2d),
  Guide.yticks(ticks=collect(-70:10:70)),
  Guide.xticks(ticks=collect(-40:10:80)))
  for i in 1:length(X)
    push!(p.layers, Gadfly.layer(x=[0.0;X[i]], y=[0.0;Y[i]], Geom.line, Gadfly.Theme(default_color=colorant"magenta"))[1])
  end
  p
end

# moved to RobotUtils.jl
# lsrBR(a) = [a[2,:];a[1,:]]';

function saveImgSeq(d::Dict{Int,Array{Float64,2}}; from::Int=1,to::Int=10,step::Int=1)
  for i in from:step:to
    p = plotLsrScanFeats(lsrBR(d[i]));
    Gadfly.draw(PNG(string("imgs/img",i,".png"),25cm,25cm),p)
  end
  nothing
end



# --------------------------------------------------------------
# transfered in from IncrementalInference

## TODO -- you were here with port starboard lines
function stbPrtLineLayers!(pl, Xpp, Ypp, Thpp; l::Float64=5.0)
    if DISABLESTBPRTLINES
      return nothing
    end


    lnstpr = [0.0;l;0.0]
    lnstpg = [0.0;-l;0.0]

    Rd  =SE2(lnstpr)
    Gr = SE2(lnstpg)

    for i in 1:length(Xpp)
      lnstt = [Xpp[i];Ypp[i];Thpp[i]]
      Ps = SE2(lnstt)
      lnr = se2vee(Ps*Rd)
      lng = se2vee(Ps*Gr)
      xsr = [Xpp[i];lnr[1]]
      ysr = [Ypp[i];lnr[2]]
      xsg = [Xpp[i];lng[1]]
      ysg = [Ypp[i];lng[2]]

      push!(pl.layers, layer(x=xsr, y=ysr, Geom.path(), Gadfly.Theme(default_color=colorant"red", line_width=1.5pt))[1] )
      push!(pl.layers, layer(x=xsg, y=ysg, Geom.path(), Gadfly.Theme(default_color=colorant"green", line_width=1.5pt))[1] )
    end
    nothing
end

# function lblsFromTo(from,to)
#   lbls=String[]
#   [push!(lbls, "$(i)") for i in from:to]
#   return lbls
# end

function drawPoses(fg::FactorGraph; from::Int=0,to::Int=99999999,
                    meanmax=:max, lbls=true, drawhist=true,
                    spscale::Float64=5.0,
                    api::DataLayerAPI=IncrementalInference.localapi  )
    #Gadfly.set_default_plot_size(20cm, 30cm)
    Xp,Yp = get2DPoseSamples(fg, from=from, to=to)
    Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; LBLS=String[];
    if meanmax == :mean
      Xpp,Ypp, Thpp, LBLS = get2DPoseMeans(fg, from=from, to=to, api=api)
    elseif meanmax == :max
      Xpp,Ypp, Thpp, LBLS = get2DPoseMax(fg, from=from, to=to, api=api)
    end

    # lbls = lblsFromTo(1,length(Xpp))
    psplt = Union{}
    if lbls
      psplt = Gadfly.plot(
      Gadfly.layer(x=Xpp,y=Ypp,label=LBLS,Geom.path(), Theme(line_width=1pt), Geom.label)
      )
    else
      psplt = Gadfly.plot(
      Gadfly.layer(x=Xpp,y=Ypp,Geom.path(), Theme(line_width=1pt))
      )
    end
    stbPrtLineLayers!(psplt, Xpp, Ypp, Thpp, l=spscale)
    if drawhist
      push!(psplt.layers,  Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)[1] )#(xbincount=100, ybincount=100))
    end
    psplt
end

function drawLandms(fg::FactorGraph;
              from::Int=0, to::Int=99999999, minnei::Int=0,
              meanmax=:max,
              lbls=true,showmm=false,drawhist=true,
              c="red",
              MM=Union{},
              api::DataLayerAPI=IncrementalInference.localapi  )
    #Gadfly.set_default_plot_size(20cm, 30cm)
    Xp,Yp = get2DLandmSamples(fg, from=from, to=to)
    Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; lblstags=String[];
    if meanmax==:mean
      Xpp,Ypp, t, lbltags = get2DLandmMeans(fg, from=from, to=to, api=api)
    elseif meanmax==:max
      Xpp,Ypp, t, lbltags = get2DLandmMax(fg, from=from, to=to,showmm=showmm,MM=MM, api=api)
    end

    if lbls
      psplt = Gadfly.plot(
      Gadfly.layer(x=Xpp,y=Ypp, label=lbltags, Geom.point, Theme(line_width=1pt, default_color=parse(Colorant,c), point_size=1pt), Geom.label)
      # ,Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)#(xbincount=100, ybincount=100)
      )
    else
      psplt = Gadfly.plot(
      Gadfly.layer(x=Xpp,y=Ypp, Geom.point, Theme(line_width=1pt, default_color=parse(Colorant,c), point_size=1pt))
      )
    end

    if drawhist
      push!(psplt.layers, Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)[1])#(xbincount=100, ybincount=100)
    end

    psplt
end

function drawPosesLandms(fg::FactorGraph;
                    from::Int=0, to::Int=99999999, minnei::Int=0,
                    meanmax=:max,lbls=true,drawhist=true, MM=Union{}, showmm=true,
                    spscale::Float64=5.0,
                    api::DataLayerAPI=IncrementalInference.localapi  )
  p = drawPoses(fg, from=from,to=to,meanmax=meanmax,lbls=lbls,drawhist=drawhist, spscale=spscale, api=api)
  pl = drawLandms(fg, from=from, to=to, minnei=minnei,lbls=lbls,drawhist=drawhist, MM=MM, showmm=showmm, api=api)
  for l in pl.layers
    push!(p.layers, l)
  end
  return p
end

function drawSubmaps(fgl::FactorGraph, fromto::Array{Int,2};
                    m1hist=false,m2hist=false,m3hist=false, showmm=false, MM=Union{},
                    api::DataLayerAPI=IncrementalInference.localapi )
  p = drawLandms(fgl, from=fromto[1,1], to=fromto[1,2], drawhist=m1hist, showmm=showmm, MM=MM, api=api)
  if size(fromto,1) >1
    p2 = drawLandms(fgl, from=fromto[2,1], to=fromto[2,2], drawhist=m2hist,c="blue", showmm=showmm, MM=MM, api=api)
    for l in p2.layers
      push!(p.layers, l)
    end
  end
  if size(fromto,1) >2
    p3 = drawLandms(fgl, from=fromto[3,1], to=fromto[3,2], drawhist=m3hist,c="magenta", showmm=showmm, MM=MM, api=api)
    for l in p3.layers
      push!(p.layers, l)
    end
  end
  return p
end

function drawSubmaps(fgl::FactorGraph, fromto::Array{Int,1}; spread::Int=25,
                    m1hist=false,m2hist=false,m3hist=false, showmm=false, MM=Union{})
  ft = zeros(Int,length(fromto),2)
  for i in 1:length(fromto)
    ft[i,1] = fromto[i]-spread; ft[i,2] = fromto[i]+spread;
  end
  drawSubmaps(fgl, ft, m1hist=m1hist, m2hist=m2hist, m3hist=m3hist, showmm=showmm, MM=MM)
end

# function getKDEMax(p::BallTreeDensity;N=200)
#   m = zeros(p.bt.dims)
#   for i in 1:p.bt.dims
#     mm = marginal(p,[i])
#     rangeV = getKDERange(mm)
#     X = linspace(rangeV[1],rangeV[2],N)
#     yV = evaluateDualTree(mm,X)
#     m[i] = X[findfirst(yV,maximum(yV))]
#   end
#   return m
# end

function investigatePoseKDE(p::BallTreeDensity, p0::BallTreeDensity)
    # co = ["black"; "blue"]
    # h = Union{}
    # x = plotKDE([marginal(p,[1]); marginal(p0,[1])], c=co )
    # y = plotKDE([marginal(p,[2]); marginal(p0,[2])], c=co )
    # if p.bt.dims >= 3
    #   th = plotKDE([marginal(p,[3]); marginal(p0,[3])], c=co )
    #   h = hstack(x,y,th)
    # else
    #   h = hstack(x,y)
    # end
    #
    # return h
    return investigateMultidimKDE(p, p0)
end


function investigatePoseKDE(p::Array{BallTreeDensity,1})
    # co = ["black"; "blue"; "green"; "red"; "magenta"; "cyan"; "cyan1"; "cyan2";
    # "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta";
    # "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"]
    # # compute all the marginals
    # Pm = Array{Array{BallTreeDensity,1},1}()
    # push!(Pm,stackMarginals(p,1)) #[marginal(p[1],[1]); marginal(p[2],[1])]
    # push!(Pm,stackMarginals(p,2)) #[marginal(p[1],[2]); marginal(p[2],[2])]
    #
    # h = Union{}
    # x = plotKDE(Pm[1], c=co )
    # y = plotKDE(Pm[2], c=co )
    # if p[1].bt.dims >= 3
    #   #Pm3 = [marginal(p[1],[3]); marginal(p[2],[3])]
    #   push!(Pm,stackMarginals(p,3)) # [marginal(p[1],[3]); marginal(p[2],[3])]
    #   th = plotKDE(Pm[3], c=co )
    #   h = hstack(x,y,th)
    # else
    #   h = hstack(x,y)
    # end
    # return h
    return investigateMultidimKDE(p)
end

function investigatePoseKDE(p::BallTreeDensity)
    # x = plotKDE(marginal(p,[1]) )
    # y = plotKDE(marginal(p,[2]) )
    # if p.bt.dims >= 3
    #   th = plotKDE(marginal(p,[3]) )
    #   return hstack(x,y,th)
    # end
    # return hstack(x,y)
    return investigateMultidimKDE(p)
end


function drawMarginalContour(fgl::FactorGraph, lbl::String;
    xmin=-150,xmax=150,ymin=-150,ymax=150,n=200,
    api::DataLayerAPI=IncrementalInference.localapi )
  #
  p = getVertKDE(fgl,lbl, api=api)  # p = getKDE(getVert(fgl,lbl))
  Gadfly.plot(z=(x,y)->evaluateDualTree(p,vectoarr2([x,y]))[1],
    x=linspace(xmin,xmax,n),
    y=linspace(ymin,ymax,n),
    Geom.contour,
    Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    Guide.title(lbl)
  )
end

function accumulateMarginalContours(fgl, order;
    xmin=-150,xmax=150,ymin=-150,ymax=150,n=200,
    api::DataLayerAPI=IncrementalInference.localapi )
  #
  pl = drawMarginalContour(fgl, order[1],xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,n=n, api=api)
  pl2 = nothing
  PL = []
  for or in order[1:end]
    pl2 = drawMarginalContour(fgl, or, xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,n=n, api=api)
    push!(PL, pl2)
    push!(pl.layers, pl2.layers[1])
  end
  return pl, PL
end




function plotPose3Pairs(fgl::FactorGraph, sym::Symbol; fill::Bool=true)
  p1= plotKDE(fgl, :x1, dims=[1;2], fill=fill)
  p2 = plotKDE(fgl, :x1, dims=[6;3], fill=fill)
  p3 = plotKDE(fgl, :x1, dims=[4;5], fill=fill)
  Gadfly.draw(PDF("/tmp/RoMEvstackPose3.pdf",15cm, 20cm), vstack(p1,p2,p3) )
  @async run(`evince /tmp/RoMEvstackPose3.pdf`)
  nothing
end













#
