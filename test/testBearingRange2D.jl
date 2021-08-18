using RoME
using Statistics
# , Distributions
using Test
using DistributedFactorGraphs
using TensorCast
import Base: convert

##

@testset "test sampling from BearingRange factor..." begin

##

p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))

fg = initfg()
addVariable!(fg, :x0, Pose2)
addVariable!(fg, :x1, Point2)
addFactor!(fg, [:x0;:x1], p2br, graphinit=false)

meas = sampleFactor(IIF._getCCW(fg, :x0x1f1), 100)

##

# meas = getSample(p2br, 100)
M = getManifold(p2br)
mcords = vee.(Ref(M), Ref(identity_element(M)), meas[1])

mu = Statistics.mean(mcords)
sigma = Statistics.std(mcords)

@test abs(mu[1]) < 0.1
@test 0.05 < abs(sigma[1]) < 0.2

@test abs(mu[2] - 20.0) < 1.0
@test 0.5 < abs(sigma[2]) < 1.5

##

end


@testset "test BearingRange factor residual function..." begin

##

# dummy variables
fg = initfg()
X0 = addVariable!(fg, :x0, Pose2)
X1 = addVariable!(fg, :x1, Point2)

##

p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))

xi = getPointIdentity(Pose2)
li = zeros(2); li[1] = 20.0;

M = getManifold(p2br)
_zi = [0,20.0]
zi = (Manifolds.hat(M, identity_element(M), _zi),)

res = calcFactorResidualTemporary( p2br, (Pose2, Point2), zi, (xi, li)) 
#
# calcFactorResidualTemporary(p2br, zi, (Pose2, xi), (Point2, li))

@show res
@test norm(res) < 1e-14

##

xi = getPointIdentity(Pose2)
li = zeros(2); li[2] = 20.0;
_zi = [pi/2,20.0]
zi = (Manifolds.hat(M, identity_element(M), _zi),)
# idx = 1
# res = zeros(2)
# p2br(res, fmd, idx, zi, xi, li)

res = calcFactorResidualTemporary( p2br, (Pose2, Point2), zi, (xi, li)) 

@show res
@test norm(res) < 1e-14

Xi = zeros(3); Xi[3] = pi/2
xi = getPoint(Pose2, Xi) 
li = zeros(2); li[2] = 20.0;
_zi = [0.0,20.0]
zi = (Manifolds.hat(M, identity_element(M), _zi),)

res = calcFactorResidualTemporary( p2br, (Pose2, Point2), zi, (xi, li))

#
@show res
@test norm(res) < 1e-14

##

Xi = zeros(3); Xi[3] = -pi/2
xi = getPoint(Pose2, Xi) 
li = zeros(2); li[1] = 20.0;
# zi = ([0.0;pi/2],[0.0;20.0],)
_zi = [pi/2,20.0]
zi = (Manifolds.hat(M, identity_element(M), _zi),)
# idx = 2
# res = zeros(2)
# p2br(res, fmd, idx, zi, xi, li)

res = calcFactorResidualTemporary( p2br, (Pose2, Point2), zi, (xi, li))


@show res
@test norm(res) < 1e-14

##
# #TODO Update to new CalcFactor 
# # test parametric Pose2Point2BearingRange
# f = Pose2Point2BearingRange(Normal(0.0,1), Normal(10.0,1))
# @test isapprox(f([0.,0,0], [10.,0]), 0, atol = 1e-9)
# @test isapprox(f([0,0,pi/2], [0.,10]), 0, atol = 1e-9)

# f = Pose2Point2BearingRange(Normal(pi/2,1), Normal(10.0,1))
# @test isapprox(f([0.,0,0], [0.,10]), 0, atol = 1e-9)
# @test isapprox(f([0,0,pi/2], [-10.,0]), 0, atol = 1e-9)

# f = Pose2Point2BearingRange(Normal(pi,1), Normal(10.0,1))
# @test isapprox(f([0.,0,0], [-10.,0]), 0, atol = 1e-9)
# @test isapprox(f([0,0,pi/2], [0.,-10]), 0, atol = 1e-9)

# f = Pose2Point2BearingRange(Normal(-pi/2,1), Normal(10.0,1))
# @test isapprox(f([0.,0,0], [0.,-10]), 0, atol = 1e-9)
# @test isapprox(f([0,0,pi/2], [10.,0]), 0, atol = 1e-9)

end




@testset "test unimodal bearing range factor, solve for landmark..." begin

##

# Start with an empty graph
N = 1
fg = initfg()

#add pose with partial constraint
addVariable!(fg, :x0, Pose2)
addFactor!(fg, [:x0], PriorPose2(MvNormal(zeros(3), 0.01*Matrix{Float64}(LinearAlgebra.I, 3,3))), graphinit=false)
# force particular initialization
setVal!(fg, :x0, [getPointIdentity(Pose2)])

##----------- sanity check that predictbelief plumbing is doing the right thing
_pts, = predictbelief(fg, :x0, ls(fg, :x0), N=75)
@cast pts[j,i] := DFG.getCoordinates.(Pose2, _pts)[i][j]
@test sum(abs.(Statistics.mean(pts,dims=2)) .< [0.1; 0.1; 0.1]) == 3
@test sum([0.05; 0.05; 0.05] .< Statistics.std(pts,dims=2) .< [0.15; 0.15; 0.15]) == 3
#------------

# Add landmark
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
li = zeros(2); li[1] = 20.0;
setVal!(fg, :l1, [li])


# Add bearing range measurement between pose and landmark
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br, graphinit=false)

# there should be just one (the bearingrange) factor connected to :l1
@test length(ls(fg, :l1)) == 1
# drawGraph(fg, show=true)

# check the forward convolution is working properly
_pts, = predictbelief(fg, :l1, ls(fg, :l1), N=75)
@cast pts[j,i] := _pts[i][j]
@test sum(abs.(Statistics.mean(pts,dims=2) - [20.0; 0.0]) .< [2.0; 2.0]) == 2
@test sum([0.1; 0.1] .< Statistics.std(pts,dims=2) .< [3.0; 3.0]) == 2

# using Gadfly, KernelDensityEstimate, KernelDensityEstimatePlotting
#
# pl = plotKDE(kde!(pts))
# pl.coord = Coord.Cartesian(xmin=-5,xmax=25, ymin=-10.0,ymax=10)
# pl

##

end


@testset "test unimodal bearing range factor, solve for pose..." begin

##

# Start with an empty graph
N = 1
fg = initfg()

# Add landmark
addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
addFactor!(fg, [:l1], PriorPoint2(MvNormal([20.0;0.0], Matrix(Diagonal([1.0;1.0].^2)))),  graphinit=false ) # could be IIF.Prior
li = zeros(2); li[1] = 20.0;
setVal!(fg, :l1, [li])

#add pose with partial constraint
addVariable!(fg, :x0, Pose2)
# force particular initialization
setVal!(fg, :x0, [getPointIdentity(Pose2)])

# Add bearing range measurement between pose and landmark
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br, graphinit=false)

# there should be just one (the bearingrange) factor connected to :l1
@test length(ls(fg, :x0)) == 1
# writeGraphPdf(fg)

# check the forward convolution is working properly
_pts, = predictbelief(fg, :x0, ls(fg, :x0), N=75)
p_μ = mean(SpecialEuclidean(2), _pts)

_pts = getCoordinates.(Pose2, _pts)
@cast pts[j,i] := _pts[i][j]

# pts[3,:] .= TU.wrapRad.(pts[3,:])
@show abs.(Statistics.mean(pts,dims=2))
@test sum(abs.(Statistics.mean(pts,dims=2)) .< [2.0; 2.0; 2.0]) == 3
@show Statistics.std(pts,dims=2)
@test sum([0.1; 2.0; 0.01] .< Statistics.std(pts,dims=2) .< [5.0; 10.0; 2.0]) == 3

##

end

@testset "Testing Pose2Point2Bearing Initialization and Packing" begin

##

p2p2b = Pose2Point2Bearing( MvNormal([0.2,0.2,0.2], [1.0 0 0;0 1 0;0 0 1]) )
packed = convert(PackedPose2Point2Bearing, p2p2b)
p2p2bTest = convert(Pose2Point2Bearing, packed)
@test p2p2b.bearing.μ == p2p2bTest.bearing.μ
@test p2p2b.bearing.Σ.mat == p2p2bTest.bearing.Σ.mat

##

end

#=


fg = LightDFG{SolverParams}(solverParams=SolverParams())

# Add the first pose :x0
x0 = addVariable!(fg, :x0, Pose2)

p = getPoint(Pose2, [10; 10; -pi])

prior = addFactor!(fg, [:x0], PriorPose2( MvNormal([0.1,0.1,0.05]), p ) )


for i in 0:5
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/3], Matrix(Diagonal([0.1;0.1;0.1].^2))))
    addFactor!(fg, [psym;nsym], pp )
end

addVariable!(fg, :l1, Point2, tags=[:LANDMARK;])
p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x0; :l1], p2br)

p2br = Pose2Point2BearingRange(Normal(0,0.1),Normal(20.0,1.0))
addFactor!(fg, [:x6; :l1], p2br)

smtasks = Task[]
solveTree!(fg; smtasks)
#add bearing range



## ================================================================================================
## Other tests
## ================================================================================================
hist = hists[3]
step = 7
fnc_ = hist[step].f
# the data before step
csmc_ = deepcopy(hist[step].csmc);
csmc_.enableLogging = false;
csmc_.logger =  SimpleLogger(Base.devnull);#NullLogger()

##

fg = initfg()

pRight = 0.9
pWrong = 0.1
pr_noise = [0.01, 0.01, 0.001]
od_noise = [0.2; 0.2; 0.2]
lm_noise = [0.01, 0.01, 0.001]

#x0 prior
addVariable!(fg, :x0, Pose2)
prpo = MvNormal([0,0,-pi], pr_noise)
addFactor!(fg, [:x0], PriorPose2(MvNormal(rand(prpo), pr_noise)))

#l1 and l2
addVariable!(fg, :l1, Pose2, tags=[:LANDMARK])
addVariable!(fg, :l2, Pose2, tags=[:LANDMARK])

#x0 to l1 or l2
p2ln = MvNormal([1, 1, pi], lm_noise)
p2p = Pose2Pose2(MvNormal(rand(p2ln), lm_noise))
addFactor!(fg, [:x0; :l1; :l2], p2p, multihypo = [1.0, pRight, pWrong])

#x0 to x1
addVariable!(fg, :x1, Pose2)
pp = MvNormal([1.0,0,0], od_noise)
addFactor!(fg, [:x0,:x1], Pose2Pose2(MvNormal(rand(pp), od_noise)))

#x1 to l1 or l2
p2ln = MvNormal([0, 1, pi], lm_noise)
p2p = Pose2Pose2(MvNormal(rand(p2ln), lm_noise))
addFactor!(fg, [:x1; :l1; :l2], p2p, multihypo = [1.0, pRight, pWrong])

#x1 to l2 or l1
p2ln = MvNormal([0, -1, pi], lm_noise)
p2p = Pose2Pose2(MvNormal(rand(p2ln), lm_noise))
addFactor!(fg, [:x1; :l2; :l1], p2p, multihypo = [1.0, pRight, pWrong])

#aslo good example without x2
#x1 to x2
addVariable!(fg, :x2, Pose2)
pp = MvNormal([1.0,0,0], od_noise)
addFactor!(fg, [:x1,:x2], Pose2Pose2(MvNormal(rand(pp), od_noise)))

#x2 to l1 or l2
p2ln = MvNormal([-1, 1, pi], lm_noise)
p2p = Pose2Pose2(MvNormal(rand(p2ln), lm_noise))
addFactor!(fg, [:x2; :l1; :l2], p2p, multihypo = [1.0, pRight, pWrong])

#x2 to l2 or l1
p2ln = MvNormal([-1, -1, pi], lm_noise)
p2p = Pose2Pose2(MvNormal(rand(p2ln), lm_noise))
addFactor!(fg, [:x2; :l2; :l1], p2p, multihypo = [1.0, pRight, pWrong])


##

fg.solverParams.inflation=0.1
fg.solverParams.spreadNH=0.1
solveTree!(fg)

plotSLAM2D(fg)




##

# Random.seed!(42) # The answer to reproducable noise

fg = LightDFG(solverParams=SolverParams(graphinit=false, gibbsIters=5, spreadNH=5.0))

pRight = 0.8
pWrong = 0.2
pr_noise = 0.01
od_noise = 0.1
lm_noise = 0.01

# true positions
# x0 at 0
# x1 at 1
# l1 at 1
# l2 at 2

#x0 prior
addVariable!(fg, :x0, ContinuousScalar)
prpo = Normal(0.0, pr_noise)
addFactor!(fg, [:x0], Prior(Normal(rand(prpo), pr_noise)))

#l1 and l2
addVariable!(fg, :l1, ContinuousScalar, tags=[:LANDMARK])
addVariable!(fg, :l2, ContinuousScalar, tags=[:LANDMARK])

#x0 to l1 or l2
p2ln = Normal(1.0, lm_noise)
p2p = LinearRelative(Normal(rand(p2ln), lm_noise))
addFactor!(fg, [:x0; :l1; :l2], p2p, multihypo = [1, pRight, pWrong])
# addFactor!(fg, [:x0; :l1], p2p) #this one used for sanity check

#x0 to x1
addVariable!(fg, :x1, ContinuousScalar)
pp = Normal(1.0, od_noise)
addFactor!(fg, [:x0,:x1], LinearRelative(Normal(rand(pp), od_noise)))

#x1 to l1 or l2
p2ln = Normal(0.0, lm_noise)
p2p = LinearRelative(Normal(rand(p2ln), lm_noise))
addFactor!(fg, [:x1; :l1; :l2], p2p, multihypo = [1, pRight, pWrong])
# addFactor!(fg, [:x1; :l1], p2p) #this one used for sanity check

#x1 to l2 or l1
p2ln = Normal(1.0, lm_noise)
p2p = LinearRelative(Normal(rand(p2ln), lm_noise))
addFactor!(fg, [:x1; :l2; :l1], p2p, multihypo = [1, pRight, pWrong])
# addFactor!(fg, [:x1; :l2], p2p) #this one used for sanity check

#

# prescribe an elimination order to get a single clique
eo = [:l2,:x1,:x0,:l1]
# fg.solverParams.graphinit=true

fg.solverParams.inflation=3.0
fg.solverParams.spreadNH=1.0
fg.solverParams.dbg = true
##
smtasks = Task[]
tree, _, = solveTree!(fg; smtasks, eliminationOrder=eo) #, smtasks=smtasks, recordcliqs=ls(fg));


# hists = fetchCliqHistoryAll!(smtasks)

plotKDE(fg, ls(fg))

##

@test isapprox(DFG.getPPESuggested(fg, :x0)[], 0, atol = 0.2) 
@test isapprox(DFG.getPPESuggested(fg, :x1)[], 1, atol = 0.2) 
@test isapprox(DFG.getPPESuggested(fg, :l1)[], 1, atol = 0.2) 

L2 = getBelief(fg, :l2)
L2_ = manikde!(2 .+ 0.1*randn(size(getPoints(L2),2)), ContinuousScalar)

# test that there is at least a mode present
mmd(L2_, L2, ContinuousScalar)
@test isapprox(DFG.getPPESuggested(fg, :l2)[], 2, atol = 0.2) 



=#
