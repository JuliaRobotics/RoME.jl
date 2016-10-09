using RoME, IncrementalInference
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


GTp = Dict{ASCIIString, Vector{Float64}}()
GTl = Dict{ASCIIString, Vector{Float64}}()


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

function landmsInRange(GTl::Dict{ASCIIString, Vector{Float64}}, cur::Vector{Float64};
			lim::Float64=100.0)
	inrange = Dict{ASCIIString, Float64}()
	for gtl in GTl
		rho = norm(gtl[2]-cur)
		if rho < lim
			inrange[gtl[1]] = rho
		end
	end
	return inrange
end

function isInFG!(fgl::FactorGraph, lbl::ASCIIString; N=100, ready=0)
	v = nothing
	if !haskey(fgl.IDs, lbl)
		init = 300*randn(2,N)
		v = addNode!(fgl, lbl, init, diagm([1000.0;1000.0]), N=N, ready=ready)
	else
		v = getVert(fgl, lbl)
	end
	v
end

function addLandmsOnPose!(fgl::FactorGraph, pose::Graphs.ExVertex, GTl::Dict{ASCIIString, Float64};
			ready=0,N=100)
	for gtl in GTl
		println("addLandmsOnPose! -- adding $(gtl[1])")
		v = isInFG!(fgl, gtl[1], N=N,ready=ready)
		# add the constraint
		ppr = Point2DPoint2DRange([gtl[2]], 2.0, [1.0])
		addFactor!(fgl, [pose;v], ppr, ready=ready)
	end
	nothing
end

function addNewPose!(fgl::FactorGraph, from::ASCIIString, lbl::ASCIIString, GTp;
			ready=0, N=N)

	init = 300*randn(2,N)

	v = addNode!(fgl, lbl, init, diagm([1000.0;1000.0]), N=N, ready=ready)
	rhoZ = norm(GTp[lbl]-GTp[from])
	ppr = Point2DPoint2DRange([rhoZ], 3.0, [1.0])
	f = addFactor!(fgl, [getVert(fgl, from);v], ppr, ready=ready)
	pts = evalFactor2(fgl, f, v.index)
	setVal!(v, pts)
	updateFullVert!(fgl,v)
	getVert(fgl, lbl)
end

function drive(fgl::FactorGraph, GTp, GTl, from, to; N=100)
	v = addNewPose!(fgl, from, to, GTp, N=N)
	# v = getVert(fgl, to)
	addLandmsOnPose!(fgl, v, landmsInRange(GTl, GTp[to], lim=120.0), N=N )
	writeGraphPdf(fgl)
	nothing
end

function batchsolve(fgl::FactorGraph)
	tree = wipeBuildNewTree!(fgl, drawpdf=true)
	inferOverTree!(fgl, tree, N=N)
	nothing
end

# function drawMarginalContour(fgl::FactorGraph, lbl::ASCIIString;xmin=-150,xmax=150,ymin=-150,ymax=150,n=200)
# 	p = getKDE(getVert(fgl,lbl))
# 	plot(z=(x,y)->evaluateDualTree(p,([x;y]')')[1],
# 		x=linspace(xmin,xmax,n),
# 		y=linspace(ymin,ymax,n),
# 		Geom.contour,
# 		Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
# 	)
# end

N = 300
fg = initfg()
fg.sessionname="SESSranges"
# fg.cg = cloudGraph

# Some starting position
init = 300*randn(2,N)
v1 = addNode!(fg, "l100", init, diagm([1000.0;1000.0]), N=N, ready=0)


lmv1 = landmsInRange(GTl, GTp["l100"])
addLandmsOnPose!(fg, v1, lmv1, N=N )

# must pin landmarks for guage
pp2 = PriorPoint2D(GTl["l1"], diagm([1.0;1.0]), [1.0])
f = addFactor!(fg,[getVert(fg,"l1")], pp2)
pp2 = PriorPoint2D(GTl["l2"], diagm([1.0;1.0]), [1.0])
f = addFactor!(fg, [getVert(fg,"l2")], pp2)



writeGraphPdf(fg)

batchsolve(fg)
drawLandms(fg,showmm=true)
drawLandms(fg,showmm=true,from=100,to=100)
drawLandms(fg,showmm=true,from=3,to=3)

RoME.drawMarginalContour(fg,"l3")

addNewPose!(fg, "l100", "l101", GTp, N=N)
writeGraphPdf(fg)

batchsolve(fg)
drawLandms(fg,showmm=true)
drawLandms(fg,showmm=true,from=101,to=101)

RoME.drawMarginalContour(fg,"l3")
RoME.drawMarginalContour(fg,"l101")



addLandmsOnPose!(fg, getVert(fg, "l101"),
								landmsInRange(GTl, GTp["l101"]), N=N )

writeGraphPdf(fg)

batchsolve(fg)
drawLandms(fg,showmm=true)

RoME.drawMarginalContour(fg,"l3")
RoME.drawMarginalContour(fg,"l100")


drive(fg, GTp, GTl, "l101", "l102", N=N)

batchsolve(fg)
drawLandms(fg,showmm=true)





drive(fg, GTp, GTl, "l102", "l103", N=N)

batchsolve(fg)
drawLandms(fg,showmm=true)






drive(fg, GTp, GTl, "l103", "l104", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)





drive(fg, GTp, GTl, "l104", "l105", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)





drive(fg, GTp, GTl, "l105", "l106", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)





drive(fg, GTp, GTl, "l106", "l107", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)



drive(fg, GTp, GTl, "l107", "l108", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)



drive(fg, GTp, GTl, "l108", "l109", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)




drive(fg, GTp, GTl, "l109", "l110", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)



drive(fg, GTp, GTl, "l110", "l111", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)



drive(fg, GTp, GTl, "l111", "l112", N=N)
batchsolve(fg)
drawLandms(fg,showmm=true)



#
