## add more julia processes
using Distributed
nprocs() < 4 ? addprocs(4-nprocs()) : nothing

# access modules/namespaces
using RoME
using RoMEPlotting

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
addVariable!(fg, :l100, Point2)

# add three landmarks
addVariable!(fg, :l1, Point2)
addVariable!(fg, :l2, Point2)
addVariable!(fg, :l3, Point2)

# and put priors on :l101 and :l102
addFactor!(fg, [:l1;], PriorPoint2(MvNormal(GTl[:l1], Matrix{Float64}(LinearAlgebra.I,2,2))) )
addFactor!(fg, [:l2;], PriorPoint2(MvNormal(GTl[:l2], Matrix{Float64}(LinearAlgebra.I,2,2))) )


# first range measurement
rhoZ1 = norm(GTl[:l1]-GTp[:l100])
ppr = Point2Point2Range( Normal(rhoZ1, 2.0) )
addFactor!(fg, [:l100;:l1], ppr)

# second range measurement
rhoZ2 = norm(GTl[:l2]-GTp[:l100])
ppr = Point2Point2Range( Normal(rhoZ2, 3.0) )
addFactor!(fg, [:l100; :l2], ppr)

# second range measurement
rhoZ3 = norm(GTl[:l3]-GTp[:l100])
ppr = Point2Point2Range( Normal(rhoZ3, 3.0) )
addFactor!(fg, [:l100; :l3], ppr)


## solve system
tree, smt, hist = solveTree!(fg)


# plot the first results
plotKDE(fg, [:l1;:l2], dims=[1;2])

plotKDE(fg, :l100, dims=[1;2], levels=6)

drawLandms(fg, from=1, to=101)# |> PDF("/tmp/test.pdf"); #@async run(`evince /tmp/test.pdf`)


function vehicle_drives_to!(fgl::G, pos_sym::Symbol, GTp::Dict, GTl::Dict; measurelimit::R=150.0) where {G <: AbstractDFG, R <: Real}
  currvar = ls(fgl)
  prev_sym = Symbol("l$(maximum(Int[parse(Int,string(currvar[i])[2:end]) for i in 2:length(currvar)]))")
  if !(pos_sym in currvar)
    println("Adding variable vertex $pos_sym, not yet in fgl::FactorGraph.")
    addVariable!(fgl, pos_sym, Point2)
    @show rho = norm(GTp[prev_sym] - GTp[pos_sym])
    ppr = Point2Point2Range( Normal(rho, 3.0) )
    addFactor!(fgl, [prev_sym;pos_sym], ppr)
  else
    @warn "Variable node $pos_sym already in the factor graph."
  end
  beacons = keys(GTl)
  for ll in beacons
    rho = norm(GTl[ll] - GTp[pos_sym])
    # Check for feasible measurements:  vehicle within 150 units from the beacons/landmarks
    if rho < measurelimit
      ppr = Point2Point2Range( Normal(rho, 3.0) )
      if !(ll in currvar)
        println("Adding variable vertex $ll, not yet in fgl::FactorGraph.")
        addVariable!(fgl, ll, Point2)
      end
      addFactor!(fgl, [pos_sym;ll], ppr)
    end
  end
  nothing
end




vehicle_drives_to!(fg, :l101, GTp, GTl)
vehicle_drives_to!(fg, :l102, GTp, GTl)


tree, smt, hist = solveTree!(fg)


pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 0:2], dims=[1;2])

pl = plotKDE(fg, [:l3;:l4], dims=[1;2], levels=4)



vehicle_drives_to!(fg, :l103, GTp, GTl)
vehicle_drives_to!(fg, :l104, GTp, GTl)


tree, smt, hist = solveTree!(fg)

pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 0:4], dims=[1;2]) |> PDF("/tmp/test.pdf")
@async run(`evince /tmp/test.pdf`)


# drive further
vehicle_drives_to!(fg, :l105, GTp, GTl)
vehicle_drives_to!(fg, :l106, GTp, GTl)

tree, smt, hist = solveTree!(fg)

pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 0:4], dims=[1;2]) |> PDF("/tmp/test.pdf")
@async run(`evince /tmp/test.pdf`)



vehicle_drives_to!(fg, :l107, GTp, GTl)

tree, smt, hist = solveTree!(fg)




vehicle_drives_to!(fg, :l108, GTp, GTl)

tree, smt, hist = solveTree!(fg)

pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 2:8], dims=[1;2], levels=6)
pl |> PDF("/tmp/test.pdf"); @async run(`evince /tmp/test.pdf`)



vehicle_drives_to!(fg, :l109, GTp, GTl)
vehicle_drives_to!(fg, :l110, GTp, GTl)

tree, smt, hist = solveTree!(fg)


vehicle_drives_to!(fg, :l111, GTp, GTl)
vehicle_drives_to!(fg, :l112, GTp, GTl)

tree, smt, hist = solveTree!(fg)


pl = plotKDE(fg, [Symbol("l$(100+i)") for i in 7:12], dims=[1;2])
pl |> PDF("/tmp/test.pdf"); @async run(`evince /tmp/test.pdf`)


pl = plotKDE(fg, [:l1;:l2;:l3;:l4], dims=[1;2], levels=4)
pl |> PDF("/tmp/test.pdf"); @async run(`evince /tmp/test.pdf`)


pl = drawLandms(fg)
pl |> PDF("/tmp/test.pdf"); @async run(`evince /tmp/test.pdf`)


drawTree(tree, show=true, imgs=true)
