using RoME
using Test
using TensorCast
using DistributedFactorGraphs
using Plots

function scatDir!(pts; arrow_scale=1.0, dir_color=:cyan, kwargs...)

    points = map(eachcol(pts)) do p
        u = cos(p[3]) * arrow_scale
        v = sin(p[3]) * arrow_scale
        x = p[1]
        y = p[2]
        plot!([x,x+u], [y,y+v], lc=dir_color, la=0.9)
        (x,y)
    end
    scatter!(points; kwargs...)
end

##

lmp_noise = diagm([0.01;0.01;0.1].^2)

# new factor graph
fg = initfg()
fg.solverParams.N = 100
fg.solverParams.inflation = 0.0
# fg.solverParams.inflation = 0.0
# landmarks
addVariable!(fg, :l1, Pose2)
addVariable!(fg, :l2, Pose2)
addVariable!(fg, :l3, Pose2)

addFactor!(fg, [:l1], PriorPose2(MvNormal([-10.0;1.0-10.0; 0],lmp_noise)), graphinit=false)
addFactor!(fg, [:l2], PriorPose2(MvNormal([-10.0+sqrt(3)/2;-0.5-10.0; 0],lmp_noise)), graphinit=false)
addFactor!(fg, [:l3], PriorPose2(MvNormal([-10.0-sqrt(3)/2;-0.5-10.0; 0],lmp_noise)), graphinit=false)

initAll!(fg)
# pose
addVariable!(fg, :x1, Point2)
addFactor!(fg, [:l1, :x1], Pose2Point2Bearing(Normal(-pi/2,0.001)), graphinit=false)
addFactor!(fg, [:l2, :x1], Pose2Point2Bearing(Normal(-pi-pi/6,0.001)), graphinit=false)
addFactor!(fg, [:l3, :x1], Pose2Point2Bearing(Normal(pi/6,0.001)), graphinit=false)


# initManual!(fg, :x1, [0.01*randn(2,100);-randn(1,100)])
# initManual!(fg, :x1, [30.0*randn(2,100);randn(1,100)])

init_points = [samplePoint(TranslationGroup(2), MvNormal([0.,0], diagm([2.,2].^2))) for i=1:fg.solverParams.N]
initManual!(fg, :x1, init_points)


##

points = getPoints(fg, :l1)
@cast pts[j,i] :=  getCoordinates.(Pose2, points)[i][j]
scatter(pts[1,:], pts[2,:], markersize=10,size=(600,600), legend=false)
scatDir!(pts; arrow_scale, markersize=1, dir_color=:red)

points = getPoints(fg, :l2)
@cast pts[j,i] := getCoordinates.(Pose2, points)[i][j]
scatter!(pts[1,:], pts[2,:], markersize=10)
scatDir!(pts; arrow_scale, markersize=1, dir_color=:blue)

points = getPoints(fg, :l3)
@cast pts[j,i] :=  getCoordinates.(Pose2, points)[i][j]
scatter!(pts[1,:], pts[2,:], markersize=10)
scatDir!(pts; arrow_scale, markersize=1, dir_color=:green)

##
arrow_scale = 0.1
fg.solverParams.inflation = 3.0
mkd, dens, labels, sum_inferred_dim = localProduct(fg, :x1)

points = getPoints(dens[1])
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=2)

points = getPoints(dens[2])
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=2)

points = getPoints(dens[3])
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=2)


# points = sample(mkd, 300)[1]
points = getPoints(mkd)
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=2)
# scatDir!(pts; arrow_scale=0.15, markersize=2, dir_color=:magenta)


##

initManual!(fg, :x1, points)


##

fg.solverParams.inflation = 1.0
eliminationOrder = [:x1, :l1, :l2, :l3]
tree = solveTree!(fg; eliminationOrder)
points = getPoints(getBelief(fg, :x1))
@cast pts[j,i] := points[i][j]
pts = collect(pts)

scatter!(pts[1,:], pts[2,:], markersize=2)
##