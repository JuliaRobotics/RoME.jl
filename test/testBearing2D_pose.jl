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

# =================================================================================


## basic point fixed pose loose
arrow_scale = 1.2
fg = initfg()
fg.solverParams.inflation = 0.0
# fg.solverParams.inflation = 1.0
# fg.solverParams.inflation = 3.0

addVariable!(fg, :x1, Pose2)
# addFactor!(fg, [:x1], PriorPose2(MvNormal([0.0; 0; 0], diagm([1;1;1*3.14].^2))))
# addFactor!(fg, [:x1], PriorPose2(MvNormal([10.0; 0; 0], diagm([1;1;1*3.14].^2))))
# addFactor!(fg, [:x1], PriorPose2(MvNormal([0.0; 0; 0], diagm([3;3;1*3.14].^2))))

grd = range(-12,12, step=1.5)
@info "N is $(length(grd)^2)"

points = ProductRepr{Tuple{Vector{Float64}, Matrix{Float64}}}[]
for i=grd,j=grd
    s,c = sincos(rand()*2pi)
    push!(points, ProductRepr([i,j], [c -s; s c]))
end

initManual!(fg, :x1, points)

addVariable!(fg, :l1, Point2)
addFactor!(fg, [:l1], PriorPoint2(MvNormal([7.07,7.07], diagm([0.01,0.01].^2))))

initAll!(fg)

points = getPoints(fg, :l1)
@cast pts[j,i] := points[i][j]
scatter(pts[1,:], pts[2,:], markersize=7, size=(600,600), legend=false)

points = getPoints(fg, :x1)
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
scatter!(pts[1,:], pts[2,:], markersize=2)
# scatter!(3*cos.(pts[3,:]), 3*sin.(pts[3,:]), aspect_ratio=:equal)


addFactor!(fg, [:x1;:l1], Pose2Point2Bearing(Normal(0, 0.01)))

points = approxConv(fg, :x1l1f1, :x1)
# initManual!(fg, :x1, points)

@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
# scatter!(pts[1,:], pts[2,:], aspect_ratio=:equal, markersize=2)
# # scatter!(5*cos.(pts[3,:]), 5*sin.(pts[3,:]), aspect_ratio=:equal)
# heads = map(eachcol(pts)) do p
#     cos(p[3])*arrow_scale, sin(p[3])*arrow_scale
# end
# points = map(eachcol(pts)) do p
#     (p[1],p[2])
# end

# quiver!(points, quiver=heads, color=:cyan)
scatDir!(pts; arrow_scale, markersize=2, dir_color=:cyan)

##
addVariable!(fg, :l2, Point2)
addFactor!(fg, [:l2], PriorPoint2(MvNormal([-7.07,-7.07], diagm([0.01,0.01].^2))))

initAll!(fg)

points = getPoints(fg, :l2)
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:],markersize=7)

addFactor!(fg, [:x1;:l2], Pose2Point2Bearing(Normal(pi/2, 0.01)))

points = approxConv(fg, :x1l2f1, :x1)
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
scatDir!(pts; arrow_scale, markersize=2, dir_color=:magenta)
# scatter!(pts[1,:], pts[2,:], aspect_ratio=:equal)
# scatter!(5*cos.(pts[3,:]), 5*sin.(pts[3,:]), aspect_ratio=:equal)
# heads = map(eachcol(pts)) do p
#     cos(p[3])*arrow_scale, sin(p[3])*arrow_scale
# end
# points = map(eachcol(pts)) do p
#     (p[1],p[2])
# end

# quiver!(points, quiver=heads, color=:magenta)

##

points = predictbelief(fg, :x1, ls(fg, :x1))[1]
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
scatDir!(pts; arrow_scale, markersize=2, dir_color=:lightgreen)

##
# Using local product

points = getPoints(fg, :l1)
@cast pts[j,i] := points[i][j]
scatter(pts[1,:], pts[2,:], markersize=10,size=(600,600), legend=false)

points = getPoints(fg, :l2)
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=10)

fg.solverParams.inflation = 0.5
mkd, dens, labels, sum_inferred_dim = localProduct(fg, :x1)


points = getPoints(dens[1])
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
# scatDir!(pts; arrow_scale, markersize=1, dir_color=:magenta)

points = getPoints(dens[2])
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
# scatDir!(pts; arrow_scale, markersize=2, dir_color=:cyan)


allpoints = vcat(getPoints.(dens)...)

mkd = manikde!(getManifold(Pose2()), vcat(getPoints.(dens)...))
points = sample(mkd, 100)[1]
# points = getPoints(mkd)
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
scatDir!(pts; arrow_scale, markersize=2, dir_color=:red)


##
p1 = nothing
for lw = [3,10,20]
θ = range(-3*pi/4, stop = pi/4, length = 100)
x = 10* cos.(θ)
y = 10* sin.(θ)
p1 = plot!(x, y, linewidth=lw, α = 0.1, color=:green)
end
p1
##
## ========================================================================================================
points = getPoints(fg, :l1)
@cast pts[j,i] := points[i][j]
scatter(pts[1,:], pts[2,:], markersize=10,size=(600,600), legend=false)

points = getPoints(fg, :l2)
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=10)

##
solveTree!(fg)
points = getPoints(fg, :x1)
@cast pts[j,i] := getCoordinates.(Pose2,points)[i][j]
scatDir!(pts; arrow_scale, markersize=2, dir_color=:cyan)

##

## ========================================================================================================
using MultiKDE
using Distributions, Random, Plots

##
dims = [ContinuousDim(), ContinuousDim(), ContinuousDim()]
bw = [1, 1, 0.1]

pts = getCoordinates.(Pose2, allpoints)

kd = KDEMulti(dims, bw, pts)

grd = range(-12,12, step=1)

##
grid = Vector{Vector{Float64}}()
ps = Float64[]
ps2d = Vector{Tuple{Float64,Float64,Float64}}()
ii = 1
jj = 1
kk = 1
angle_range = -pi:0.1:pi
data_3d = zeros(length(grd),length(grd),length(angle_range))
for (ii,i)=enumerate(grd), (jj,j)=enumerate(grd), (kk,k)=enumerate(angle_range)
    push!(grid,[i,j,k])
    p = MultiKDE.pdf(kd, [i,j,k])
    push!(ps, p)
    push!(ps2d, (i,j,p))

    data_3d[ii,jj,kk] = p


end

##

idxs = findall(x->x>2.1, data_3d)
coords = zeros(3,length(idxs))
for (i,idx) in enumerate(idxs)
    a = grd[idx[1]]
    b = grd[idx[2]]
    c = angle_range[idx[3]]
    coords[:,i] = [a,b,c]
end

##
points = getPoints(fg, :l1)
@cast pts[j,i] := points[i][j]
scatter(pts[1,:], pts[2,:], markersize=10,size=(600,600), legend=false)

points = getPoints(fg, :l2)
@cast pts[j,i] := points[i][j]
scatter!(pts[1,:], pts[2,:], markersize=10)

scatDir!(coords; arrow_scale, markersize=2, dir_color=:green)


##






