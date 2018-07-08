## add more julia processes
# nprocs() < 7 ? addprocs(7-nprocs()) : nothing

# tell Julia that you want to use these modules/namespaces
using RoME, Distributions


GTp = Dict{Symbol, Vector{Float64}}()
GTp[:l100] = [0.0;0]
GTp[:l101] = [50.0;0]
GTp[:l102] = [100.0;0]
GTp[:l103] = [100.0;50.0]
GTp[:l104] = [100.0;100.0]
GTp[:l105] = [50.0;100.0]
GTp[:l106] = [0.0;100.0]
GTp[:l107] = [0.0;50.0]
GTp[:l108] = [0.0;-50.0]
GTp[:l109] = [0.0;-100.0]
GTp[:l110] = [50.0;-100.0]
GTp[:l111] = [100.0;-100.0]
GTp[:l112] = [100.0;-50.0]

GTl = Dict{Symbol, Vector{Float64}}()
GTl[:l1] = [10.0;30]
GTl[:l2] = [30.0;-30]
GTl[:l3] = [80.0;40]
GTl[:l4] = [120.0;-50]



# create the factor graph object
fg = initfg()

# first pose with no initial estimate
addNode!(fg, :l100, Point2)

# add three landmarks
addNode!(fg, :l1, Point2)
addNode!(fg, :l2, Point2)
addNode!(fg, :l3, Point2)

# and put priors on :l101 and :l102
addFactor!(fg, [:l1;], PriorPoint2D(GTl[:l1], eye(2), [1.0]))
addFactor!(fg, [:l2;], PriorPoint2D(GTl[:l2], eye(2), [1.0]))


# first range measurement
rhoZ1 = norm(GTl[:l1]-GTp[:l100])
ppr = Point2DPoint2DRange([rhoZ1], 2.0, [1.0])
addFactor!(fg, [:l100;:l1], ppr)

# second range measurement
rhoZ2 = norm(GTl[:l2]-GTp[:l100])
ppr = Point2DPoint2DRange([rhoZ2], 3.0, [1.0])
addFactor!(fg, [:l100; :l2], ppr)

# second range measurement
rhoZ3 = norm(GTl[:l3]-GTp[:l100])
ppr = Point2DPoint2DRange([rhoZ3], 3.0, [1.0])
addFactor!(fg, [:l100; :l3], ppr)


writeGraphPdf(fg)


tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



using RoMEPlotting


drawLandms(fg, from=100)

drawLandms(fg, to=99)



using KernelDensityEstimatePlotting, Gadfly



pl = plotKDE(fg, :l100, dims=[1;2])
Gadfly.draw(PDF("/tmp/testL100.pdf", 20cm, 10cm),pl)

pl = plotKDE(fg, [:l1;:l2], dims=[1;2], levels=4)
Gadfly.draw(PDF("/tmp/testL1_2.pdf", 20cm, 10cm),pl)





0


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
# init = 300*randn(2,N)
v1 = addNode!(fg, :l1, Point2, N=N, ready=0)

# Two landmarks
L1, L2, L3 = [10.0;30], [30.0;-30], [70.0;30]
l1 = addNode!(fg, :l2, Point2, N=N, ready=0)
l2 = addNode!(fg, :l3, Point2, N=N, ready=0)

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
# v2 = addOdoFG!(fg, :l4, [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P1-P2)], 3.0, [1.0])
v2 = addNode!(fg, :l4, Point2, N=N, ready=0)
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
l3 = addNode!(fg, :l5, Point2, N=N, ready=0)
addFactor!(fg, [v2;l3], ppr, ready=0)





#pose 3 in truth is at
P3 = [100.0;0.0]

# drive forward 50 units
# v2 = addOdoFG!(fg, :l4, [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P2-P3)], 3.0, [1.0])
v3 = addNode!(fg, :l6, Point2, N=N, ready=0)
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
# v2 = addOdoFG!(fg, :l4, [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P3-P4)], 3.0, [1.0])
v4 = addNode!(fg, :l7, Point2, N=N, ready=0)
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
# v2 = addOdoFG!(fg, :l4, [50.0;0.0;0.0], odoCov, ready=0, N=N)
ppr = Point2DPoint2DRange([norm(P4-P5)], 3.0, [1.0])
v5 = addNode!(fg, :l8, Point2, N=N, ready=0)
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
