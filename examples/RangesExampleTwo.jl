# video at: https://vimeo.com/190052649


addprocs(7)
using RoME, IncrementalInference, Gadfly, Colors, KernelDensityEstimate
# for p in procs()
# # @everywhere begin
# remotecall_fetch(p, ()->using RoME)
# remotecall_fetch(p, ()->using IncrementalInference)
# remotecall_fetch(p, ()->using Gadfly)
# remotecall_fetch(p, ()->using Colors)
# # end
# end

@everywhere using RoME, RoMEPlotting
@everywhere using IncrementalInference
@everywhere using Gadfly
@everywhere using Colors
@everywhere using KernelDensityEstimate
# end

# , Gadfly

# using CloudGraphs, Neo4j
# include("BlandAuthDB.jl")
#
# configuration = CloudGraphs.CloudGraphConfiguration(dbaddress, 7474, dbusr, dbpwd, mongoaddress, 27017, false, "", "");
# cloudGraph = connect(configuration);
# # register types of interest in CloudGraphs
# registerGeneralVariableTypes!(cloudGraph)
# IncrementalInference.setCloudDataLayerAPI!()

# togglePrtStbLines()


GTp = Dict{String, Vector{Float64}}()
GTl = Dict{String, Vector{Float64}}()


GTp["l100"] = [0.0;0]
GTp["l101"] = [50.0;0]
GTp["l102"] = [100.0;0]
GTp["l103"] = [100.0;50.0]
GTp["l104"] = [100.0;100.0]
GTp["l105"] = [50.0;100.0]
GTp["l106"] = [0.0;100.0]
GTp["l107"] = [0.0;50.0]
GTp["l108"] = [0.0;-50.0]
GTp["l109"] = [0.0;-100.0]
GTp["l110"] = [50.0;-100.0]
GTp["l111"] = [100.0;-100.0]
GTp["l112"] = [100.0;-50.0]

GTl["l1"] = [10.0;30]
GTl["l2"] = [30.0;-30]
GTl["l3"] = [80.0;40]
GTl["l4"] = [120.0;-50]

cursor = [0.0;0]

function landmsInRange(GTl::Dict{String, Vector{Float64}}, cur::Vector{Float64};
			lim::Float64=100.0)
	inrange = Dict{String, Float64}()
	for gtl in GTl
		rho = norm(gtl[2]-cur)
		if rho < lim
			inrange[gtl[1]] = rho
		end
	end
	return inrange
end

function isInFG!(fgl::FactorGraph, lbl::Symbol; N=100, ready=0)
	v = nothing
	if !haskey(fgl.IDs, lbl)
		init = 300*randn(2,N)
		v = addNode!(fgl, lbl, init, diagm([1000.0;1000.0]), N=N, ready=ready)
	else
		v = getVert(fgl, lbl)
	end
	v
end

function addLandmsOnPose!(fgl::FactorGraph, pose::Graphs.ExVertex, GTl::Dict{String, Float64};
			ready=0,N=100)
	for gtl in GTl
		println("addLandmsOnPose! -- adding $(gtl[1])")
		v = isInFG!(fgl, Symbol(gtl[1]), N=N,ready=ready)
		# add the constraint
		ppr = Point2DPoint2DRange([gtl[2]], 2.0, [1.0])
		addFactor!(fgl, [pose;v], ppr, ready=ready)
	end
	nothing
end

function addNewPose!(fgl::FactorGraph, from::Symbol, lbl::Symbol, GTp;
			ready=0, N=N)

	init = 300*randn(2,N)

	v = addNode!(fgl, lbl, init, diagm([1000.0;1000.0]), N=N, ready=ready)
	rhoZ = norm(GTp[string(lbl)]-GTp[string(from)])
	ppr = Point2DPoint2DRange([rhoZ], 3.0, [1.0])
	f = addFactor!(fgl, [from,lbl], ppr, ready=ready)
  # f = addFactor!(fgl, [getVert(fgl, from);v], ppr, ready=ready)
  initializeNode!(fgl,lbl)
	# pts = evalFactor2(fgl, f, v.index)
	# setVal!(v, pts)
	# updateFullVert!(fgl,v)
	getVert(fgl, lbl)
end

function drive(fgl::FactorGraph, GTp, GTl, from, to; N=100)
	v = addNewPose!(fgl, from, to, GTp, N=N)
	# v = getVert(fgl, to)
	addLandmsOnPose!(fgl, v, landmsInRange(GTl, GTp[string(to)], lim=120.0), N=N )
  println("added landmark")
	writeGraphPdf(fgl)
	nothing
end

function batchsolve(fgl::FactorGraph; N::Int=100)
	tree = wipeBuildNewTree!(fgl, drawpdf=true)
	inferOverTree!(fgl, tree, N=N)
	nothing
end

@everywhere begin
function drawQuadLandms(fgl; file="", w=30cm, h=20cm,
			h11 = Gadfly.context(),
			h12 = Gadfly.context())
	h21 = haskey(fgl.IDs,"l3") ? drawMarginalContour(fgl,"l3") : Gadfly.context()
	h22 = haskey(fgl.IDs,"l4") ? drawMarginalContour(fgl,"l4") : Gadfly.context()
	hh = vstack(hstack(h11,h12),hstack(h21,h22))
	if length(file) > 0
		Gadfly.draw(PNG(file,w,h),hh)
	end
	hh
end


function drawLandmMargOver(fgl::FactorGraph, lbl::String,
		file,
		gtlayers;
		w=30cm, h=20cm)

	hover = Gadfly.context()
	if haskey(fgl.IDs,lbl)
		hover =  drawMarginalContour(fgl,lbl)
		# for gtl in gtlayers
		# 	push!(hover.layers, gtl)
		# end
	end

	if length(file) > 0
		Gadfly.draw(PNG(file,w,h), hover)
	end
	hover
end


function drawGroundTruth(GTp, orderp, GTl=nothing, orderl=[]; drawranges=true, interp=false, t=0)
	LAYERS = Gadfly.Layer[]
	X = Float64[]
	Y = Float64[]
	for o in orderp
		push!(X, GTp[o][1])
		push!(Y, GTp[o][2])
	end
	if interp
		X[end] = t*GTp[orderp[length(orderp)-1]][1]+(1-t)*GTp[orderp[end]][1]
		Y[end] = t*GTp[orderp[length(orderp)-1]][2]+(1-t)*GTp[orderp[end]][2]
	end
	push!(LAYERS, Gadfly.layer(x=X, y=Y, Geom.point)[1])
	for i in 2:(length(orderp))
		push!(LAYERS, Gadfly.layer(x=[X[i-1];X[i]], y=[Y[i-1];Y[i]], Geom.path())[1])
	end
	Xlr = Float64[]
	Ylr = Float64[]
	Lblr = String[]
	Xla = Float64[]
	Yla = Float64[]
	Lbla = String[]
	for o in orderl
		if o == "l1" || o == "l2"
			push!(Lblr, o)
			push!(Xlr, GTl[o][1])
			push!(Ylr, GTl[o][2])
		else
			push!(Lbla, o)
			push!(Xla, GTl[o][1])
			push!(Yla, GTl[o][2])
		end
	end
	if length(orderl) > 0
		push!(LAYERS, Gadfly.layer(x=Xlr, y=Ylr, Geom.point, Theme(default_color=colorant"red", point_size=4pt)  )[1]) # label=Lblr, Geom.label
		push!(LAYERS, Gadfly.layer(x=Xla, y=Yla, Geom.point, Theme(default_color=colorant"cyan", point_size=4pt) )[1]) # label=Lbla, Geom.label
		if drawranges
			idx = !interp ? length(X) : length(X)-1
			gtlcur = landmsInRange(GTl, [X[idx];Y[idx]])
			for lbl in keys(gtlcur)
				push!(LAYERS, Gadfly.layer(x=[X[idx];GTl[lbl][1]],y=[Y[idx];GTl[lbl][2]],Geom.path(), Theme(default_color=colorant"magenta"))[1])
			end
		end
	end
	Gadfly.plot(LAYERS,
		Guide.title("Ground Truth"),
		Coord.Cartesian(xmin=-150,xmax=150,ymin=-150, ymax=150))
end

function drawConvenience1(fgl, ii, fr, folderloc, h12, h12a, interp)
	drawQuadLandms(fgl,file="$(folderloc)/overlms$(fr).png",h12=h12);
	draw(PNG("$(folderloc)/gt$(fr).png",30cm,20cm),h12)
	ii -= interp ? 1 : 0
	lstPosePl = RoME.drawMarginalContour(fgl,"l$(100+ii)")
	ii += interp ? 1 : 0
	draw(PNG("$(folderloc)/lst$(fr).png",30cm,20cm),lstPosePl)
	push!(h12a.layers, lstPosePl.layers[1])
	draw(PNG("$(folderloc)/lstOver$(fr).png",30cm,20cm),h12a)
	nothing
end

function drawConvenience2(fgl, fr, folderloc, h12b)
	ppl = drawLandms(fgl,showmm=true,lbls=false,from=100)
	ppl.coord = Coord.Cartesian(xmin=-150,xmax=150,ymin=-150,ymax=150)

	@sync begin
		allposesfile2 = "$(folderloc)/allPoses$(fr).png"
		@async draw(PNG(allposesfile2,30cm,20cm),ppl)

		# for l in ppl.layers
		templyrs = union(h12b.layers, ppl.layers)
		h12b.layers = templyrs
		# end

		allposesfile = "$(folderloc)/allPosesOver$(fr).png"
		draw(PNG(allposesfile,30cm,20cm),h12b)
	end
	nothing
end

end #everywhere

function drawAllVidImages(GTp, GTl, fgl, ii, fr; drawranges=true, interp=false, t=0)
	poselbls = ["l$(100+i)" for i in 0:ii]
	lmlbls = ["l$(i)" for i in 1:4]
	folderloc = "/home/dehann/irosVid/new"
	rr = Dict{Int,Future}()

	h12 = drawGroundTruth(GTp, poselbls, GTl, lmlbls, drawranges=drawranges, interp=interp, t=t)
	# overlay ground truth and landmark estimation
	rr[1] = @spawn drawLandmMargOver(fgl, "l3", "$(folderloc)/lm3Over$(fr).png", h12.layers)
	rr[2] = @spawn drawLandmMargOver(fgl, "l4", "$(folderloc)/lm4Over$(fr).png", h12.layers)

	h12a = drawGroundTruth(GTp, poselbls, GTl, lmlbls, drawranges=drawranges, interp=interp, t=t)
	rr[3] = @spawn drawConvenience1(fgl, ii, fr, folderloc, h12, h12a, interp)

	h12b = drawGroundTruth(GTp, poselbls, GTl, lmlbls, drawranges=drawranges, interp=interp, t=t)
	rr[4] = @spawn drawConvenience2(fgl, fr, folderloc, h12b)

	for r in rr
		fetch(r[2])
	end

	# adding kld plots
  nothing
end


function draw30AllFast(GTp, GTl, fgl, th, offset, start, ending;
				interp=true, drawranges=true)
	# @sync begin
	# 	@async [drawAllVidImages(GTp, GTl, fgl, th, i+offset, interp=interp, drawranges=drawranges, t=(30-i)/30.0) for i in start:minimum([10;ending])];
	# 	@async [drawAllVidImages(GTp, GTl, fgl, th, i+offset, interp=interp, drawranges=drawranges, t=(30-i)/30.0) for i in 11:minimum([20;ending])];
	# 	@async [drawAllVidImages(GTp, GTl, fgl, th, i+offset, interp=interp, drawranges=drawranges, t=(30-i)/30.0) for i in 21:ending];
	# end

  # doesnt work properly yet, redrawing old
	@sync begin
		for i in start:ending
			@async drawAllVidImages(GTp, GTl, fgl, th, i+offset, interp=interp, drawranges=drawranges, t=(30-i)/30.0)
		end
	end

	nothing
end


function layerCircle(;cent=[0.0;0.0], radius=1.0, c="deepskyblue", N=200)
	TH = linspace(0,2pi,200)
	X = real(radius*exp(TH*im))+cent[1]
	Y = imag(radius*exp(TH*im))+cent[2]
	layer(x=X,y=Y, Geom.path(), Theme(default_color=parse(Colorant,c)))
end
function plotCircle(;cent=[0.0;0.0], radius=1.0, c="deepskyblue", N=200, drawcent=false)
	PL = []
	push!(PL, layerCircle(cent=cent, radius=radius, c=c, N=N) )
	!drawcent ? nothing : push!(PL, layer(x=[cent[1]],y=[cent[2]], Geom.point, Theme(default_color=parse(Colorant,c))))
	Gadfly.plot(PL...)
end


function drawFirstPoseIllustration(GTl, GTp)
	l1, l2 = GTl["l1"], GTl["l2"]
	p1, p2 = GTp["l100"], GTp["l101"]
	rho1, rho2 = norm(p1-l1), norm(p1-l2)

	PL = Gadfly.Layer[]

	PL = union(PL, plotCircle(cent=[0.0;0],radius=3.0, c="blue", drawcent=false).layers)
	PL = union(PL, plotCircle(cent=[36.0;12],radius=3.0, c="blue", drawcent=false).layers)

	PL = union(PL, plotCircle(cent=l1,radius=rho1, c="red", drawcent=true).layers)
	PL = union(PL, plotCircle(cent=l2,radius=rho2, c="red", drawcent=true).layers)

	Gadfly.plot(PL,
		Guide.title("Parametric illustration"),
		Coord.Cartesian(xmin=-150,xmax=150,ymin=-150, ymax=150) )
end


function drawSecondPoseIllustration(GTl, GTp)
	l1, l2 = GTl["l1"], GTl["l2"]
	p1, p2 = GTp["l100"], GTp["l101"]
	rho1, rho2 = norm(p1-l1), norm(p1-l2)

	PL = Gadfly.Layer[]

	PL = union(PL, plotCircle(cent=[0.0;0],radius=norm(p1-p2), c="blue", drawcent=true).layers)
	PL = union(PL, plotCircle(cent=[36.0;12],radius=norm(p1-p2), c="blue", drawcent=true).layers)


	PL = union(PL, plotCircle(cent=l1,radius=rho1, c="gray").layers)
	PL = union(PL, plotCircle(cent=l2,radius=rho2, c="gray").layers)

	Gadfly.plot(PL,
		Guide.title("Parametric illustration"),
		Coord.Cartesian(xmin=-150,xmax=150,ymin=-150, ymax=150) )
end

function drawFirstL3Illustration(GTl, GTp)
	l1, l2, l3 = GTl["l1"], GTl["l2"], GTl["l3"]
	p1, p2 = GTp["l100"], GTp["l101"]
	rho1, rho2 = norm(p1-l1), norm(p1-l2)

	PL = Gadfly.Layer[]

	PL = union(PL, plotCircle(cent=[0.0;0],radius=norm(p1-l3), c="cyan", drawcent=true).layers)
	PL = union(PL, plotCircle(cent=[36.0;12],radius=norm(p1-l3), c="cyan", drawcent=true).layers)


	PL = union(PL, plotCircle(cent=l1,radius=rho1, c="gray").layers)
	PL = union(PL, plotCircle(cent=l2,radius=rho2, c="gray").layers)

	push!(PL, layer(x=[l3[1]], y=[l3[2]], Geom.point, Theme(default_color=colorant"cyan", point_size=4pt)  )[1] )

	Gadfly.plot(PL,
		Guide.title("Parametric illustration"),
		Coord.Cartesian(xmin=-150,xmax=150,ymin=-150, ymax=150) )
end


function drawIllustrations(GTl, GTp, folderloc)
	pl = drawGroundTruth(GTp, String[], GTl, ["l$(i)" for i in 1:4], drawranges=false, interp=false, t=0)
	filename = "$(folderloc)/gtPos0.png"
	draw(PNG(filename,30cm,20cm),pl)

	# first pose location
	pl2 = drawFirstPoseIllustration(GTl, GTp)
	filename2 = "$(folderloc)/gtArc0.png"
	draw(PNG(filename2,30cm,20cm),pl2)

	pl3 = drawSecondPoseIllustration(GTl, GTp)
	filename3 = "$(folderloc)/gtPose2Pred.png"
	draw(PNG(filename3,30cm,20cm),pl3)

	pl4 = drawFirstL3Illustration(GTl, GTp)
	filename4 = "$(folderloc)/gtLm3Pred.png"
	draw(PNG(filename4,30cm,20cm),pl4)


	nothing
end


function evaluateAccuracy(fgl::FactorGraph, GTp)
  # get last pose
  xx,ll = ls(fg)
  maxPosition = 0
  for l in ll
    val = parse(Int,string(l)[2:end])
    maxPosition = val >= 100 && val > maxPosition ? val : maxPosition
  end
  vsym = Symbol("l$(maxPosition)")
  pX = getVertKDE(fg, vsym)
  truePos = GTp[string(vsym)]

  postlikeli = evaluateDualTree(pX, reshape(truePos,2,1))[1]
  MaxPointError = truePos - getKDEMax(pX, N=1000)
  return postlikeli, norm(MaxPointError)
end

function evaluateAllAccuracy(fgl::FactorGraph, GTp, GTl)
  xx,ll = ls(fg)
  LK = Float64[]
  PE = Float64[]
  for l in ll
    pX = getVertKDE(fg, l)
    truePos = haskey(GTp, string(l)) ? GTp[string(l)] : GTl[string(l)]
    push!(LK, evaluateDualTree(pX, reshape(truePos,2,1))[1])
    MaxPointError = truePos - getKDEMax(pX, N=1000)
    push!(PE, norm(MaxPointError))
  end
  return LK, PE
end


folderloc2 = "/home/dehann/irosVid/new"
drawIllustrations(GTl, GTp, folderloc2)


# N = 100

fg = initfg()

for N in collect(25:25:25)

fg = initfg()
# fg.sessionname="SESSranges"
# fg.cg = cloudGraph

# Some starting position
init = 300*randn(2,N)
v1 = addNode!(fg, :l100, init, diagm([1000.0;1000.0]), N=N, ready=0)


lmv1 = landmsInRange(GTl, GTp["l100"])
addLandmsOnPose!(fg, v1, lmv1, N=N )

# must pin landmarks for guage
pp2 = PriorPoint2D(GTl["l1"], diagm([1.0;1.0]), [1.0])
f = addFactor!(fg,[:l1], pp2)
# f = addFactor!(fg,[getVert(fg,:l1)], pp2)
pp2 = PriorPoint2D(GTl["l2"], diagm([1.0;1.0]), [1.0])
f = addFactor!(fg, [:l2], pp2)


writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","w")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","w")
write(fid, "$(pe)\n")
close(fid)


drawAllVidImages(GTp, GTl, fg, 0, 0);
draw30AllFast(GTp, GTl, fg, 1, 0, 1, 29)
# @sync begin
# @async [drawAllVidImages(GTp, GTl, fg, 1, i+0, interp=true, t=(30-i)/30.0) for i in 1:10];
# @async [drawAllVidImages(GTp, GTl, fg, 1, i+0, interp=true, t=(30-i)/30.0) for i in 11:20];
# @async [drawAllVidImages(GTp, GTl, fg, 1, i+0, interp=true, t=(30-i)/30.0) for i in 21:29];
# end
## drawAllVidImages(GTp, GTl, fg, 0, 29, interp=false)




addNewPose!(fg, :l100, :l101, GTp, N=N)
initializeNode!(fg,:l101)
writeGraphPdf(fg)
batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fg, 1, 30, 0, 14, drawranges=false, interp=false)
# [drawAllVidImages(GTp, GTl, fg, 1, 30+i, drawranges=false) for i in 0:14]






addLandmsOnPose!(fg, getVert(fg, :l101),
								landmsInRange(GTl, GTp["l101"]), N=N )
writeGraphPdf(fg)
batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fg, 1, 45, 0, 14, interp=false)
# [drawAllVidImages(GTp, GTl, fg, 1, 45+i) for i in 0:14]






fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l101, :l102, N=N)
initializeNode!(fg,:l102)
writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 2, 60, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 2, i+60, interp=true, t=(30-i)/30.0) for i in 0:29]




# initializeNode!(fg,:l102)


fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l102, :l103, N=N)
batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)
# drawLandms(fg,showmm=true)
# h12 = drawGroundTruth(GTp, [:l$(00+i)" for i in 0:3], GTl, [:l$()" for i in 1:4])
# drawQuadLandms(fg,file="irosVid/test4.png",h12=h12);
# draw(PNG("irosVid/gt4.png",30cm,20cm),h12)
# lstPosePl = RoME.drawMarginalContour(fg,:l103)
# draw(PNG("irosVid/lst4.png",30cm,20cm),lstPosePl)
# push!(h12.layers, lstPosePl.layers[1])
# draw(PNG("irosVid/lstOver4.png",30cm,20cm),h12)
# drawAllVidImages(GTp, GTl, fg, 3, 4)

draw30AllFast(GTp, GTl, fgd, 3, 90, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 3, i+90, interp=true, t=(30-i)/30.0) for i in 0:29]







fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l103, :l104, N=N)
# addLandmsOnPose!(fg, getVert(fg, :l104),
# 								landmsInRange(GTl, GTp["l104"]), N=N )
# writeGraphPdf(fg)
batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 4, 120, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 4, i+120, interp=true, t=(30-i)/30.0) for i in 0:29]







fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l104, :l105, N=N)
# addLandmsOnPose!(fg, getVert(fg, :l104),
# 								landmsInRange(GTl, GTp["l104"]), N=N )
writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 5, 150, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 5, i+150, interp=true, t=(30-i)/30.0) for i in 0:29]






fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l105, :l106, N=N)
batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

# drawAllVidImages(GTp, GTl, fgd, 6, 7)
draw30AllFast(GTp, GTl, fgd, 6, 180, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 6, i+180, interp=true, t=(30-i)/30.0) for i in 0:29]








fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l106, :l107, N=N)
batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 7, 210, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 7, i+210, interp=true, t=(30-i)/30.0) for i in 0:29]








fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l107, :l108, N=N)
# addLandmsOnPose!(fg, getVert(fg, :l108),
# 								landmsInRange(GTl, GTp["l108"]), N=N )
writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 8, 240, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 8, i+240, interp=true, t=(30-i)/30.0) for i in 0:29]







fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l108, :l109, N=N)
# addLandmsOnPose!(fg, getVert(fg, :l109),
# 								landmsInRange(GTl, GTp["l109"]), N=N )
writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 9, 270, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 9, i+270, interp=true, t=(30-i)/30.0) for i in 0:29]







fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l109, :l110, N=N)
# addLandmsOnPose!(fg, getVert(fg, :l110),
# 								landmsInRange(GTl, GTp["l110"]), N=N )
writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 10, 300, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 10, i+300, interp=true, t=(30-i)/30.0) for i in 0:29]







fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l110, :l111, N=N)
# addLandmsOnPose!(fg, getVert(fg, :l111),
# 								landmsInRange(GTl, GTp["l111"]), N=N )
writeGraphPdf(fg)

batchsolve(fg, N=N)
lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

draw30AllFast(GTp, GTl, fgd, 11, 330, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 11, i+330, interp=true, t=(30-i)/30.0) for i in 0:29]








fgd = deepcopy(fg)
drive(fg, GTp, GTl, :l111, :l112, N=N)

tsec = @elapsed batchsolve(fg, N=N)
fid = open("/home/dehann/Videos/slamedonutBatchTime_$(N).txt","w")
write(fid, "$(tsec)")
close(fid)

lk, pe = evaluateAccuracy(fg, GTp)
fid = open("/home/dehann/Videos/slamedonutLstPoseLks_$(N).txt","a")
write(fid, "$(lk)\n")
close(fid)
fid = open("/home/dehann/Videos/slamedonutLstPoseErr_$(N).txt","a")
write(fid, "$(pe)\n")
close(fid)

LK, PE = evaluateAllAccuracy(fg, GTp, GTl)

fid = open("/home/dehann/Videos/slamedonutAllPE_$(N).txt","w")
for pe in PE
  write(fid, "$(pe)\n")
end
close(fid)
fid = open("/home/dehann/Videos/slamedonutAllLK_$(N).txt","w")
for lk in LK
  write(fid, "$(lk)\n")
end
close(fid)

writeGraphPdf(fg)
wipeBuildNewTree!(fg, drawpdf=true)

draw30AllFast(GTp, GTl, fgd, 12, 360, 0, 29)
# [drawAllVidImages(GTp, GTl, fg, 12, i+360, interp=true, t=(30-i)/30.0) for i in 0:29]




draw30AllFast(GTp, GTl, fg, 12, 390, 0, 29, interp=false)
# [drawAllVidImages(GTp, GTl, fg, 12, i+390) for i in 0:29]



# run(`convert -delay 100 irosVid/test*.png irosVid/anim.gif`)

run(`ffmpeg -y -i $(folderloc2)/lstOver%d.png -threads 8 -vcodec libx264 -s 1920x1080 -b:v 2M /home/dehann/Videos/slamedonutLastPose_$(N).mp4`)
# run(`ffmpeg -y -i $(folderloc2)/lm3Over%d.png -vcodec libx264 -s 1920x1080 -b:v 2M /home/dehann/Videos/slamedonutLm3_$(N).mp4`)
# run(`ffmpeg -y -i $(folderloc2)/lm4Over%d.png -vcodec libx264 -s 1920x1080 -b:v 2M /home/dehann/Videos/slamedonutLm4_$(N).mp4`)


end



#
