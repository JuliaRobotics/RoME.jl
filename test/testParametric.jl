using RoME
using Test
using Manifolds
# using RoMEPlotting


## list of tested factors
# - PriorPose2
# - PriorPoint2
# - Pose2Pose2
# - Pose2Point2BearingRange
# - Pose2Point2
# - DynPose2VelocityPrior
# - VelPose2VelPose2

@testset "Test PriorPose2 and Pose2Pose2" begin

fg = GraphsDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))

# Add the first pose :x0
x0 = addVariable!(fg, :x0, Pose2)
prior = addFactor!(fg, [:x0], PriorPose2( MvNormal([10; 10; -pi+1e-5], [0.1;0.1;0.05])))

for i in 0:3
    psym = Symbol("x$i")
    nsym = Symbol("x$(i+1)")
    addVariable!(fg, nsym, Pose2)
    pp = Pose2Pose2(MvNormal([10.0;0;pi/2], [0.1;0.1;0.1]))
    addFactor!(fg, [psym;nsym], pp)
end
initAll!(fg)
IIF.initParametricFrom!(fg)

vars = ls(fg, r"x") 
lands = ls(fg, r"l") 
results = IIF.autoinitParametric!.(fg, [vars;lands])

# pl = plotSLAM2D(fg, solveKey=:parametric, drawContour=false, xmin=-20, xmax=20, ymin=-20, ymax=20)


PM, varLabels, r, Σ = IIF.solveGraphParametric!(fg) #autodiff=:finite)

#TODO test +-pi used pi+1e-5
#FIXME look if something is wrong with angle bounds [-pi,pi), test failed
M = getManifold(Pose2)
# @test isapprox(vardict[:x0].val, [10, 10, -pi], atol = 1e-3)
# @test isapprox(vardict[:x4].val, [10, 10, -pi], atol = 1e-3)
ϵ = getPointIdentity(M)
@test isapprox(M, r[1], exp(M, ϵ, hat(M,ϵ,[10, 10, -pi])), atol = 1e-3)
@test isapprox(M, r[2], exp(M, ϵ, hat(M,ϵ,[0, 10, -pi/2])), atol = 1e-3)
@test isapprox(M, r[3], exp(M, ϵ, hat(M,ϵ,[0, 0, 0])), atol = 1e-3)
@test isapprox(M, r[4], exp(M, ϵ, hat(M,ϵ,[10, 0, pi/2])), atol = 1e-3)
@test isapprox(M, r[5], exp(M, ϵ, hat(M,ϵ,[10, 10, -pi])), atol = 1e-3)

# IIF.updateParametricSolution(fg, vardict)
# pl = plotSLAM2D(fg; lbls=true, solveKey=:parametric, point_size=4pt, drawPoints=false, drawContour=false)
end




##
# using Random
# Random.seed!(42);
# fg = LocalDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))

# fg.solverParams.graphinit = false
# pr_noise = [0.01, 0.01, 0.001]
# od_noise = [0.2;0.2;0.2]

# σ_bearing = 0.008
# σ_range = 0.01

# addVariable!(fg, :x0, Pose2)
# prpo = MvNormal([0,0,-pi], pr_noise)
# addFactor!(fg, [:x0], PriorPose2(MvNormal(rand(prpo), pr_noise)))

# addVariable!(fg, :l1, Point2, tags=[:LANDMARK])
# addVariable!(fg, :l2, Point2, tags=[:LANDMARK])
# # p2br = Pose2Point2BearingRange(Normal(pi/4,0.01),Normal(sqrt(2),0.1))
# p2br = Pose2Point2BearingRange(Normal(pi/4 + rand(Normal(0,σ_bearing)), σ_bearing),
#                                Normal(sqrt(2) + rand(Normal(0,σ_range)), σ_range))
# addFactor!(fg, [:x0; :l1], p2br)

# addVariable!(fg, :x1, Pose2)
# pp = MvNormal([1.0,0,0], od_noise)
# addFactor!(fg, [:x0,:x1], Pose2Pose2(MvNormal(rand(pp), od_noise)))

# p2br = Pose2Point2BearingRange(Normal(pi/2 + rand(Normal(0,σ_bearing)), σ_bearing),
#                                Normal(1 + rand(Normal(0,σ_range)), σ_range))
# addFactor!(fg, [:x1; :l1], p2br)

# p2br = Pose2Point2BearingRange(Normal(-pi/2 + rand(Normal(0,σ_bearing)), σ_bearing),
#                                Normal(1 + rand(Normal(0,σ_range)), σ_range))
# addFactor!(fg, [:x1; :l2], p2br)

# addVariable!(fg, :x2, Pose2)
# pp = MvNormal([1.0,0,0], od_noise)
# addFactor!(fg, [:x1,:x2], Pose2Pose2(MvNormal(rand(pp), od_noise)))

# p2br = Pose2Point2BearingRange(Normal(3pi/4 + rand(Normal(0,σ_bearing)), σ_bearing),
#                                Normal(sqrt(2) + rand(Normal(0,σ_range)), σ_range))
# addFactor!(fg, [:x2; :l1], p2br)

# p2br = Pose2Point2BearingRange(Normal(-3pi/4 + rand(Normal(0,σ_bearing)), σ_bearing),
#                                Normal(sqrt(2) + rand(Normal(0,σ_range)), σ_range))
# addFactor!(fg, [:x2; :l2], p2br)

# initAll!(fg)

# IIF.initParametricFrom!(fg)

# vardict, result, varIds, Σ = IIF.solveGraphParametric(fg) #autodiff=:finite)
# # vardict, result, varIds, Σ = IIF.solveGraphParametric(fg, autodiff=:finite)

# IIF.updateParametricSolution(fg, vardict)

# pl = plotSLAM2D(fg; lbls=true, solveKey=:parametric, point_size=4pt, drawPoints=false, drawContour=false)




##

@testset "Test Parametric PriorPose2 and Pose2Point2" begin

fg = LocalDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x1, Pose2)
addVariable!(fg, :l1, Point2)

addFactor!(fg, [:x1], PriorPose2(MvNormal([0.,0, 0], [0.01, 0.01, 0.01])))

addFactor!(fg, [:x1; :l1], Pose2Point2(MvNormal([1.0, 1], [0.1,0.1])))

initAll!(fg)

IIF.initParametricFrom!(fg)

PM, varLabels, r, Σ = IIF.solveGraphParametric(fg) #autodiff=:finite)

M = getManifold(Pose2)
ϵ = getPointIdentity(M)

@test isapprox(M, r[1], exp(M, ϵ, hat(M,ϵ,[0, 0, 0])), atol = 1e-3)
@test isapprox(r[2], [1,  1], atol = 1e-3)

# IIF.updateParametricSolution(fg, vardict)
# pl = plotSLAM2D(fg; lbls=true, solveKey=:parametric, point_size=4pt, drawPoints=false, drawContour=false)
# push!(pl, Coord.cartesian(fixed=true),style(background_color=RGB(1,1,1)))
end


##
@testset "Test Parametric PriorPoint2 and Pose2Point2BearingRange" begin
fg = GraphsDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x1, Pose2)
addVariable!(fg, :l1, Point2)
addVariable!(fg, :l2, Point2)

addFactor!(fg, [:l1], PriorPoint2(MvNormal([1., 1], [0.01, 0.01])))
addFactor!(fg, [:l2], PriorPoint2(MvNormal([1.,-1], [0.01, 0.01])))

addFactor!(fg, [:x1; :l1], Pose2Point2BearingRange(Normal(pi/4, 0.01), Normal(sqrt(2), 0.1)))
addFactor!(fg, [:x1; :l2], Pose2Point2BearingRange(Normal(3pi/4, 0.01), Normal(sqrt(2), 0.1)))

initAll!(fg)

IIF.initParametricFrom!(fg)

PM, varLabels, r, Σ = IIF.solveGraphParametric(fg) #autodiff=:finite)

@test isapprox(SpecialEuclidean(2; vectors=HybridTangentRepresentation()), r[1], ArrayPartition([2, 0.], [0 -1; 1 0.]), atol = 1e-3)

@test isapprox(r[2], [1,  1], atol = 1e-3)
@test isapprox(r[3], [1, -1], atol = 1e-3)

# IIF.updateParametricSolution(fg, vardict)
# pl = plotSLAM2D(fg; lbls=true, solveKey=:parametric, point_size=4pt, drawPoints=false, drawContour=false)
end


@testset "Test Parametric DynPose2VelocityPrior and VelPose2VelPose2" begin

@test_broken begin
@warn "Parametric VelPose2 is broken and tests skipped"    
fg = LocalDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))

# add first pose locations
addVariable!(fg, :x0, DynPose2; nanosecondtime=0)

# Prior factor as boundary condition
pp0 = DynPose2VelocityPrior(MvNormal(zeros(3), [0.01; 0.01; 0.001]), MvNormal([10.0;0], [0.1; 0.1]))
addFactor!(fg, [:x0;], pp0)

addVariable!(fg, :x1, DynPose2;  nanosecondtime=1000_000_000)

# conditional likelihood between Dynamic Point2
dp2dp2 = VelPose2VelPose2(MvNormal([10.0;0;0], [0.01;0.01;0.001]), MvNormal([0.0;0], [0.1; 0.1]))
addFactor!(fg, [:x0;:x1], dp2dp2)

initAll!(fg)

PM, varLabels, r, Σ = IIF.solveGraphParametric!(fg)

@test isapprox(r[1], [0, 0, 0, 10, 0], atol = 1e-3)
@test isapprox(r[2], [10, 0, 0, 10, 0], atol = 1e-3)
end

end


@testset "Test Parametric PriorPoint2 and Point2Point2Range" begin

fg = LocalDFG( solverParams=SolverParams(algorithms=[:default, :parametric]))

addVariable!(fg, :x1, Point2)
addVariable!(fg, :l1, Point2)
addVariable!(fg, :l2, Point2)
addVariable!(fg, :l3, Point2)

addFactor!(fg, [:l1], PriorPoint2(MvNormal([0., 0], [0.01, 0.01])))
addFactor!(fg, [:l2], PriorPoint2(MvNormal([1., 0], [0.01, 0.01])))
addFactor!(fg, [:l3], PriorPoint2(MvNormal([0., 1], [0.01, 0.01])))

addFactor!(fg, [:x1; :l1], Point2Point2Range(Normal(sqrt(2), 0.1)))
addFactor!(fg, [:x1; :l2], Point2Point2Range(Normal(1.0, 0.1)))
addFactor!(fg, [:x1; :l3], Point2Point2Range(Normal(1.0, 0.1)))

PM, varLabels, r, Σ = IIF.solveGraphParametric(fg)

@test isapprox(r[1], [1, 1], atol = 1e-3)

end