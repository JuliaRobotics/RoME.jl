## add more julia processes
nprocs() < 4 ? addprocs(4-nprocs()) : nothing

# access modules/namespaces
using RoME, Distributions
using RoMEPlotting
using KernelDensityEstimatePlotting, Gadfly


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
addFactor!(fg, [:l1;], PriorPoint2(MvNormal(GTl[:l1], Matrix{Float64}(LinearAlgebra.I, 2,2))))
addFactor!(fg, [:l2;], PriorPoint2(MvNormal(GTl[:l2], Matrix{Float64}(LinearAlgebra.I, 2,2))))


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




# drawLandms(fg, from=100)
# drawLandms(fg, to=99)





pl = plotKDE(fg, :l100, dims=[1;2])
Gadfly.draw(PNG("/tmp/testL100.png", 20cm, 10cm),pl)

pl = plotKDE(fg, [:l1;:l2], dims=[1;2], levels=4)
Gadfly.draw(PNG("/tmp/testL1_2.png", 20cm, 10cm),pl)










function vehicle_drives_to!(fgl::FactorGraph, pos_sym::Symbol, GTp::Dict, GTl::Dict; measurelimit::R=150.0) where {R <: Real}
  currvar = union(ls(fgl)...)
  prev_sym = Symbol("l$(maximum(Int[parse(Int,string(currvar[i])[2:end]) for i in 2:length(currvar)]))")
  if !(pos_sym in currvar)
    println("Adding variable vertex $pos_sym, not yet in fgl::FactorGraph.")
    addNode!(fgl, pos_sym, Point2)
    @show rho = norm(GTp[prev_sym] - GTp[pos_sym])
    ppr = Point2DPoint2DRange([rho], 3.0, [1.0])
    addFactor!(fgl, [prev_sym;pos_sym], ppr)
  else
    @warn "Variable node $pos_sym already in the factor graph."
  end
  beacons = keys(GTl)
  for ll in beacons
    rho = norm(GTl[ll] - GTp[pos_sym])
    # Check for feasible measurements:  vehicle within 150 units from the beacons/landmarks
    if rho < measurelimit
      ppr = Point2DPoint2DRange([rho], 3.0, [1.0])
      if !(ll in currvar)
        println("Adding variable vertex $ll, not yet in fgl::FactorGraph.")
        addNode!(fgl, ll, Point2)
      end
      addFactor!(fgl, [pos_sym;ll], ppr)
    end
  end
  nothing
end






vehicle_drives_to!(fg, :l101, GTp, GTl)

# tree = wipeBuildNewTree!(fg)
# inferOverTree!(fg, tree)


vehicle_drives_to!(fg, :l102, GTp, GTl)

tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)




pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 0:2], dims=[1;2], levels=6)
Gadfly.draw(PNG("/tmp/testL100_102.png", 20cm, 10cm),pl)

pl = plotKDE(fg, [:l1;:l2], dims=[1;2], levels=4)
Gadfly.draw(PNG("/tmp/testL1_2.png", 20cm, 10cm),pl)

pl = plotKDE(fg, [:l3;:l4], dims=[1;2], levels=4)
Gadfly.draw(PNG("/tmp/testL3_4.png", 20cm, 10cm),pl)







vehicle_drives_to!(fg, :l103, GTp, GTl)

# tree = wipeBuildNewTree!(fg)
# inferOverTree!(fg, tree)


vehicle_drives_to!(fg, :l104, GTp, GTl)

tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)




pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 0:4], dims=[1;2])
Gadfly.draw(PNG("/tmp/testL100_105.png", 20cm, 10cm),pl)




vehicle_drives_to!(fg, :l105, GTp, GTl)
vehicle_drives_to!(fg, :l106, GTp, GTl)


tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)


pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 0:6], dims=[1;2], levels=6)
Gadfly.draw(PNG("/tmp/testL100_106.png", 20cm, 10cm),pl)




vehicle_drives_to!(fg, :l107, GTp, GTl)

tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



vehicle_drives_to!(fg, :l108, GTp, GTl)


tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 2:8], dims=[1;2], levels=6)
Gadfly.draw(PNG("/tmp/testL103_108.png", 20cm, 10cm),pl)






vehicle_drives_to!(fg, :l109, GTp, GTl)
vehicle_drives_to!(fg, :l110, GTp, GTl)


tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 6:10], dims=[1;2])
Gadfly.draw(PNG("/tmp/testL106_110.png", 20cm, 10cm),pl)








vehicle_drives_to!(fg, :l111, GTp, GTl)
vehicle_drives_to!(fg, :l112, GTp, GTl)


tree = wipeBuildNewTree!(fg)
inferOverTree!(fg, tree)



pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 7:12], dims=[1;2])
Gadfly.draw(PNG("/tmp/testL106_112.png", 20cm, 10cm),pl)

pl = plotKDE(fg, [:l1;:l2;:l3;:l4], dims=[1;2], levels=4)
Gadfly.draw(PNG("/tmp/testL1234.png", 20cm, 10cm),pl)



pl = drawLandms(fg, from=100)
Gadfly.draw(PNG("/tmp/testLocsAll.png", 20cm, 10cm),pl)


pl = drawLandms(fg)
Gadfly.draw(PNG("/tmp/testAll.png", 20cm, 10cm),pl)


# for ll in [:l1;:l2;:l3;:l4; :l100;:l101;:l102;:l103; :l104;:l105;:l106;:l107; :l108;:l109;:l110;:l111; :l112]
#   setVal!(fg, ll, 200*randn(2,1000))
# end


# togglePrtStbLines()

#
