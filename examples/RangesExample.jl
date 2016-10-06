using RoME, IncrementalInference

# using CloudGraphs, Neo4j
# include("BlandAuthDB.jl")
#
# configuration = CloudGraphs.CloudGraphConfiguration(dbaddress, 7474, dbusr, dbpwd, mongoaddress, 27017, false, "", "");
# cloudGraph = connect(configuration);
# # register types of interest in CloudGraphs
# registerGeneralVariableTypes!(cloudGraph)
# IncrementalInference.setCloudDataLayerAPI!()


N = 300
fg = initfg()
fg.sessionname="SESSranges"
# fg.cg = cloudGraph


# initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.1])

# Some starting position
init = 300*randn(2,N)
v1 = addNode!(fg, "l1", init, diagm([1000.0;1000.0]), N=N, ready=0)

# Two landmarks
L1, L2 = [10.0;30], [40.0;-30]
l1 = addNode!(fg, "l2", (L1')', diagm([1.0;1.0]), N=N, ready=0)
l2 = addNode!(fg, "l3", (L2')', diagm([1.0;1.0]), N=N, ready=0)

# must pin landmarks for guage
pp2 = PriorPoint2D(L1, diagm([1.0;1.0]), [1.0])
f = addFactor!(fg,[l1], pp2)
pp2 = PriorPoint2D(L2, diagm([1.0;1.0]), [1.0])
f = addFactor!(fg, [l2], pp2)

# and constraints to pose x1
rhoZ1 = norm(L1)
ppr = Point2DPoint2DRange([rhoZ1], 2.0, [1.0])
addFactor!(fg, [v1;l1], ppr, ready=0)

rhoZ2 = norm(L2)
ppr = Point2DPoint2DRange([rhoZ2], 3.0, [1.0])
addFactor!(fg, [v1;l2], ppr, ready=0)

tree = prepBatchTree!(fg, drawpdf=true)


inferOverTree!(fg, tree, N=N)
# togglePrtStbLines()

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=1,to=1)

drawLandms(fg,showmm=true,from=4)



# setDBAllReady!(fg)

# drive forward 50 units
# v2 = addOdoFG!(fg, "l4", [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([40.0], 3.0, [1.0])
v2 = addNode!(fg, "l4", init, diagm([3.0;3.0]), N=N, ready=0)
addFactor!(fg, [v1;v2], ppr, ready=0)



# get one more range
rhoZ3 = norm(L2-[40;0.0])
ppr = Point2DPoint2DRange([rhoZ3], 3.0, [1.0])
addFactor!(fg, [v2;l2], ppr, ready=0)

# solve again
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=4)
# drawPosesLandms(fg,from=2,to=2)

# setDBAllReady!(fg)


# get one more range
rhoZ4 = norm(L1-[40;0.0])
ppr = Point2DPoint2DRange([rhoZ4], 3.0, [1.0])
addFactor!(fg, [v2;l1], ppr, ready=0)


# solve again
tree = wipeBuildNewTree!(fg, drawpdf=true)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)


#
