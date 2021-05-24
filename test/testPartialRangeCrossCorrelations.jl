
# see correct result with slanted (narrow) L1 modes at y=+-5.  L1 should not be round blobs
# https://github.com/JuliaRobotics/RoME.jl/pull/434#issuecomment-817246038

using RoME
using Test


##


@testset "Test correlation induced by partials" begin

##

# start with an empty factor graph object
N = 200
fg = initfg()
getSolverParams(fg).N = N

getSolverParams(fg).useMsgLikelihoods = true

addVariable!(fg, :x0, Pose2)
addVariable!(fg, :l1, Point2)
addFactor!(fg, [:x0], PriorPose2( MvNormal([0; 0; 0], diagm([0.3;0.3;0.3].^2)) ))
ppr = Pose2Point2Range(MvNormal([7.3], diagm([0.3].^2)))
addFactor!(fg, [:x0; :l1], ppr )

addVariable!(fg, :x1, Pose2)
pp = Pose2Pose2(MvNormal([9.8;0;0.8], diagm([0.3;0.3;0.05].^2)))
ppr = Pose2Point2Range(MvNormal([6.78], diagm([0.3].^2)))
addFactor!(fg, [:x0; :x1], pp )
addFactor!(fg, [:x1; :l1], ppr )


##

tree, smt, hist = solveTree!(fg)


## check that stuff is where it should be

L1_ = getPoints(getBelief(fg, :l1))

@test 0.3*N < sum(L1_[2,:] .< 0)
@test 0.3*N < sum(0 .< L1_[2,:])

L1_n = L1_[:, L1_[2,:] .< 0]
L1_p = L1_[:, 0 .< L1_[2,:]]

mvn = fit(MvNormal, L1_n)
mvp = fit(MvNormal, L1_p)

# check diagonal structure for correlation
@test isapprox(mvn.Σ.mat[1,1], 1.7, atol=1.2)
@test isapprox(mvn.Σ.mat[2,2], 1.7, atol=1.2)
@test isapprox(mvn.Σ.mat[1,2], 1.1, atol=1.2)

@test isapprox(mvp.Σ.mat[1,1], 1.7, atol=1.2)
@test isapprox(mvp.Σ.mat[2,2], 1.7, atol=1.2)
@test isapprox(mvp.Σ.mat[1,2], -1.1, atol=1.2)

# sanity check for symmetry
@test mvn.Σ.mat - mvn.Σ.mat' |> norm < 0.01

# test means in the right location
@test isapprox(mvn.μ[1], 5.4, atol=1.0)
@test isapprox(mvn.μ[2], -4.8, atol=0.75)

@test isapprox(mvp.μ[1], 5.4, atol=1.0)
@test isapprox(mvp.μ[2], 4.8, atol=0.75)

##



end


##

# using RoMEPlotting
# Gadfly.set_default_plot_size(35cm,25cm)

##


# plotSLAM2D(fg)

# plotKDE(fg, :l1, levels=3)

##
