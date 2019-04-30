
# using Revise

##

# for drawing Bayes tree with images (see debug tricks below)
# using Cairo, Fontconfig
# using Gadfly

using RoME

#  Do some plotting
# using RoMEPlotting


## start with an empty factor graph object
fg = initfg()
fg.isfixedlag = true
fg.qfl = 15


# Add the first pose :x0
addVariable!(fg, :x0, Pose2)

# Add at a fixed location PriorPose2 to pin :x0 to a starting location (10,10, pi/4)
addFactor!(fg, [:x0], PriorPose2( MvNormal([0.0; 0.0; 0.0],
                                           Matrix(Diagonal([0.1;0.1;0.05].^2))) ), autoinit=false )

# Drive around in a hexagon
for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp, autoinit=false )
end

# Add landmarks with Bearing range measurements
addVariable!(fg, :l1, Point2, labels=["LANDMARK"])
p2br = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x0; :l1], p2br, autoinit=false )


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x6; :l1], p2br2, autoinit=false )



tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)





# Drive around in a hexagon
for i in 6:11
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp, autoinit=false )
end



# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x12; :l1], p2br2, autoinit=false )




tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)




## drive a little more

# Drive around in a hexagon
for i in 12:17
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp, autoinit=false )
end


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x18; :l1], p2br2, autoinit=false )




tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)





# Drive around in a hexagon
for i in 18:23
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp, autoinit=false )
end



# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x24; :l1], p2br2, autoinit=false )






tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)







## drive a little more

# Drive around in a hexagon
for i in 24:29
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp, autoinit=false )
end


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x30; :l1], p2br2, autoinit=false )




tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)










# Drive around in a hexagon
for i in 30:35
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;-pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp, autoinit=false )
end


# Add landmarks with Bearing range measurements
p2br2 = Pose2Point2BearingRange(Normal(0,0.03),Normal(20.0,0.5))
addFactor!(fg, [:x36; :l1], p2br2, autoinit=false )




tree = batchSolve!(fg, treeinit=true, drawpdf=true, show=true)





## Plotting



using RoMEPlotting


pl = drawPosesLandms(fg, meanmax=:max)
Gadfly.draw(Gadfly.SVG("/tmp/test2.svg"),pl)  # or PNG(...)

@async run(`eog /tmp/test2.svg`)



0



## Manually organized a non-async initialization sequence

# docliqs = [
# 14;
# 13;
# 12;
# 11;
# 4;
# 9;
# 3;
# 10;
# 8;
# 7;
# 6; # with async
# 5; # with async
# 2;
# 1
# ]
#
#
#
# for i in docliqs[1:5]
#   # alltasks[i] = @async begin
# cliq = tree.cliques[i]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true )
# if !(clst in [:upsolved; :downsolved; :marginalized])
#   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
# end
#   # end
# end # for
# # end
#
#
# drawTree(tree)
#
#
# cliq = tree.cliques[9]
# cliq.attributes["label"]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, incremental=true, limititers=1 )
#
#
#
#
# cliq = whichCliq(tree, :x1)
# cliq.attributes["label"]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, incremental=true, limititers=1 )
#
#
#
# cliq = tree.cliques[6]
# cliq.attributes["label"]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, incremental=true, limititers=1 )
#
#
# cliq = tree.cliques[7]
# cliq.attributes["label"]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, incremental=true, limititers=1 )
#
#
# cliq = tree.cliques[8]
# cliq.attributes["label"]
# clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, incremental=true, limititers=1 )










#
