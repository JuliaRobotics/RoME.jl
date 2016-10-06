using RoME, IncrementalInference

# using CloudGraphs, Neo4j
# include("BlandAuthDB.jl")
#
# configuration = CloudGraphs.CloudGraphConfiguration(dbaddress, 7474, dbusr, dbpwd, mongoaddress, 27017, false, "", "");
# cloudGraph = connect(configuration);
# # register types of interest in CloudGraphs
# registerGeneralVariableTypes!(cloudGraph)
# IncrementalInference.setCloudDataLayerAPI!()

# togglePrtStbLines()

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
L1, L2, L3 = [10.0;30], [30.0;-30], [70.0;30]
l1 = addNode!(fg, "l2", (L1')', diagm([1.0;1.0]), N=N, ready=0)
l2 = addNode!(fg, "l3", (L2')', diagm([1.0;1.0]), N=N, ready=0)

# must pin landmarks for guage
pp2 = PriorPoint2D(L1, diagm([1.0;1.0]), [1.0])
f = addFactor!(fg,[l1], pp2)
pp2 = PriorPoint2D(L2, diagm([1.0;1.0]), [1.0])
f = addFactor!(fg, [l2], pp2)

# and constraints to pose x1
P1 = [0.0;0]
rhoZ1 = norm(L1-P1)
ppr = Point2DPoint2DRange([rhoZ1], 2.0, [1.0])
addFactor!(fg, [v1;l1], ppr, ready=0)

rhoZ2 = norm(L2-P1)
ppr = Point2DPoint2DRange([rhoZ2], 3.0, [1.0])
addFactor!(fg, [v1;l2], ppr, ready=0)


#solve
writeGraphPdf(fg)
tree = prepBatchTree!(fg, drawpdf=true)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=1,to=1)
drawLandms(fg,showmm=true,from=4)


# setDBAllReady!(fg)

#pose 2 in truth is at
P2 = [50.0;0.0]

# drive forward 50 units
# v2 = addOdoFG!(fg, "l4", [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P1-P2)], 3.0, [1.0])
v2 = addNode!(fg, "l4", init, diagm([3.0;3.0]), N=N, ready=0)
addFactor!(fg, [v1;v2], ppr, ready=0)


# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=4)



# get one range to second position
rhoZ3 = norm(L2-P2)
ppr = Point2DPoint2DRange([rhoZ3], 3.0, [1.0])
addFactor!(fg, [v2;l2], ppr, ready=0)

# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=4)



# setDBAllReady!(fg)


# get another range
rhoZ4 = norm(L1-P2)
ppr = Point2DPoint2DRange([rhoZ4], 3.0, [1.0])
addFactor!(fg, [v2;l1], ppr, ready=0)


# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg, drawpdf=true)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)




# START seeing a new range signal
ppr = Point2DPoint2DRange([norm(P2-L3)], 3.0, [1.0])
l3 = addNode!(fg, "l5", init, diagm([3.0;3.0]), N=N, ready=0)
addFactor!(fg, [v2;l3], ppr, ready=0)





#pose 3 in truth is at
P3 = [100.0;0.0]

# drive forward 50 units
# v2 = addOdoFG!(fg, "l4", [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P2-P3)], 3.0, [1.0])
v3 = addNode!(fg, "l6", init, diagm([3.0;3.0]), N=N, ready=0)
addFactor!(fg, [v2;v3], ppr, ready=0)


# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=5)






# get one range to third position
rhoZ5 = norm(L2-P3)
ppr = Point2DPoint2DRange([rhoZ5], 3.0, [1.0])
addFactor!(fg, [v3;l2], ppr, ready=0)

# solve again
# writeGraphPdf(fg)
# tree = wipeBuildNewTree!(fg)

# inferOverTree!(fg, tree, N=N)
#
# drawLandms(fg,showmm=true)
#
# drawLandms(fg,showmm=true,from=5)








rhoZ6 = norm(L1-P3)
ppr = Point2DPoint2DRange([rhoZ6], 3.0, [1.0])
addFactor!(fg, [v3;l1], ppr, ready=0)





rhoZ6b = norm(L3-P3)
ppr = Point2DPoint2DRange([rhoZ6b], 3.0, [1.0])
addFactor!(fg, [v3;l3], ppr, ready=0)




# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=5,to=5)











#pose 4 in truth is at
P4 = [100.0;70.0]

# drive forward 50 units
# v2 = addOdoFG!(fg, "l4", [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P3-P4)], 3.0, [1.0])
v4 = addNode!(fg, "l7", init, diagm([3.0;3.0]), N=N, ready=0)
addFactor!(fg, [v3;v4], ppr, ready=0)


# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=5)









# get one range to fourth position
rhoZ7 = norm(L2-P4)
ppr = Point2DPoint2DRange([rhoZ7], 3.0, [1.0])
addFactor!(fg, [v4;l2], ppr, ready=0)

# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=4)







# get one range to fourth position
rhoZ8 = norm(L1-P4)
ppr = Point2DPoint2DRange([rhoZ8], 3.0, [1.0])
addFactor!(fg, [v4;l1], ppr, ready=0)

# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=6)







# get one range to fourth position
rhoZ9 = norm(L3-P4)
ppr = Point2DPoint2DRange([rhoZ9], 3.0, [1.0])
addFactor!(fg, [v4;l3], ppr, ready=0)

# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)

drawLandms(fg,showmm=true,from=5,to=5)










#pose 5 in truth is at
P5 = [150.0;70.0]

# drive forward 50 units
# v2 = addOdoFG!(fg, "l4", [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P4-P5)], 3.0, [1.0])
v5 = addNode!(fg, "l8", init, diagm([3.0;3.0]), N=N, ready=0)
addFactor!(fg, [v4;v5], ppr, ready=0)


rhoZ10 = norm(L3-P5)
ppr = Point2DPoint2DRange([rhoZ10], 3.0, [1.0])
addFactor!(fg, [v5;l3], ppr, ready=0)



rhoZ11 = norm(L2-P5)
ppr = Point2DPoint2DRange([rhoZ11], 3.0, [1.0])
addFactor!(fg, [v5;l2], ppr, ready=0)


# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)




rhoZ12 = norm(L1-P5)
ppr = Point2DPoint2DRange([rhoZ12], 3.0, [1.0])
addFactor!(fg, [v5;l1], ppr, ready=0)


# solve again
writeGraphPdf(fg)
tree = wipeBuildNewTree!(fg)

inferOverTree!(fg, tree, N=N)

drawLandms(fg,showmm=true)



#
