
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




# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=false)
# at = initInferTreeUp!(fg, tree, drawtree=true)
#
#
# ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
# # downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=true);
# downMsgPassingRecursive(ett,N=100, dbg=false, drawpdf=true);




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




tree = wipeBuildNewTree!(fg, drawpdf=true, show=true, imgs=false)



## Manually organized a non-async initialization sequence

docliqs = [
14;
13;
12;
11;
4;
9;
3;
10;
8;
7;
6; # with async
5; # with async
2;
1
]





# resetTreeCliquesForUpSolve!(tree)

# queue all the tasks
# alltasks = Vector{Task}(undef, length(treel.cliques))
# @sync begin
# i = docliqs[8]

for i in docliqs[1:7]
  # alltasks[i] = @async begin
cliq = tree.cliques[i]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true )
if !(clst in [:upsolved; :downsolved; :marginalized])
  error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
end
  # end
end # for
# end



alltasks = Vector{Task}(undef, length(tree.cliques))
# @sync begin
for i in docliqs[8:12]
  alltasks[i] = @async begin
    cliq = tree.cliques[i]
    clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=10 )
    # if !(clst in [:upsolved; :downsolved; :marginalized])
    #   error("Clique $(cliq.index), initInferTreeUp! -- cliqInitSolveUp! did not arrive at the desired solution statu: $clst")
    # end
  end
end # for
# end
#
# alltasks


i = docliqs[13]
cliq = tree.cliques[i]
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true,  limititers=1)
0



cliq = tree.cliques[6]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )


cliq = tree.cliques[5]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )




cliq = tree.cliques[2]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )



cliq = tree.cliques[1]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )





# fg_bu = deepcopy(fg)
# tree_bu = deepcopy(tree)
# drawTree(tree_bu)
# drawTree(tree)


ts7 = @async begin

cliq = tree.cliques[7]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

end


ts8 = @async begin

cliq = tree.cliques[8]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )

end


# why are the arrays empty??
getData(cliq).upInitMsgs



cliq = tree.cliques[10]
cliq.attributes["label"]
getCliqStatus(cliq)
clst = cliqInitSolveUp!(fg, tree, cliq, drawtree=true, limititers=1 )












0









at = initInferTreeUp!(fg, tree, drawtree=true, limititers=100)




ett = ExploreTreeType(fg, tree, tree.cliques[1], nothing, NBPMessage[])
# downMsgPassingIterative!(ett,N=100, dbg=false, drawpdf=true);
downMsgPassingRecursive(ett,N=100, dbg=false, drawpdf=true);



0




using RoMEPlotting


pl = drawPosesLandms(fg, meanmax=:max)
Gadfly.draw(Gadfly.SVG("/tmp/test2.svg"),pl)  # or PNG(...)
@async run(`eog /tmp/test2.svg`)



#
