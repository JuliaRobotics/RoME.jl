#FIXME update to manifolds
#FIXME add to tests
# test simulated PPE in reverse


using Test
using RoME
using DistributedFactorGraphs
using Manifolds: hat, exp
##

@testset "ensure solveParametricBinary is working" begin

##


fg = initfg()

addVariable!(fg, :x0, Pose2)
pp = PriorPose2( MvNormal(zeros(3), diagm(0.01*ones(3)) ))
addFactor!(fg, [:x0;], pp)


# reference ppe on :x0
refVal = zeros(3)
refKey = :simulated
ppe = DFG.MeanMaxPPE(refKey, refVal, refVal, refVal)
setPPE!(fg[:x0], refKey, DFG.MeanMaxPPE, ppe)

pp2 = Pose2Pose2( MvNormal([0;0;-pi+0.01], diagm(0.03*ones(3)) ))

##

M = getManifold(Pose2())
ϵ = getPointIdentity(Pose2())

X = hat(M, ϵ, [0;0;-pi]) #measurement
p = ϵ # variable from
q = ϵ  # variable to
@test isapprox( abs.(calcFactorResidualTemporary(pp2, (Pose2, Pose2), X, (p,  q))), [0;0;pi] )

q = exp(M, ϵ, hat(M, ϵ, [0;0;-pi])) # variable to
@test isapprox( calcFactorResidualTemporary(pp2, (Pose2, Pose2), X, (p,  q)), [0;0;0], atol=1e-14 )

q = exp(M, ϵ, hat(M, ϵ, [0;0;pi])) # variable to
@test isapprox( calcFactorResidualTemporary(pp2, (Pose2, Pose2), X, (p,  q)), [0;0;0], atol=1e-14 )


# wXjhat = SE2(zeros(3))*SE2([0;0;-pi])
# jXjhat = SE2(wxj) \ wXjhat
# return se2vee(jXjhat)


## do in factor graph

addVariable!(fg, :x1, Pose2)
addFactor!( fg, [:x0; :x1], pp2, graphinit=false )


##

refVal = accumulateFactorMeans(fg, [:x0f1; :x0x1f1])

@test isapprox(refVal[1:2], [0;0], atol=1e-4)
@test 0.9pi < abs(refVal[3])


## alternate API

isAlready, simPPE, genLabel = IIF._checkVariableByReference(fg, :x0, r"x\\d+", Pose2, pp2)

@test isapprox(simPPE.suggested[1:2], [0;0], atol=1e-2)
@test 0.9pi < abs(simPPE.suggested[3])


##

end


@testset "Secondary binary parametric solve test case" begin

##

fg = initfg()

addVariable!(fg, :x2, Pose2)
addFactor!(fg, [:x2], 
            PriorPose2( MvNormal( [
                        15.000000000016204;
                        8.660254037814505;
                        2.0943951023931953; #+1e-5;
                      ], 
                      diagm([0.01;0.01;0.01]) )) )
#

addVariable!(fg, :x3, Pose2)

pp = Pose2Pose2(MvNormal([10;0;pi/3], diagm([0.01;0.01;0.01])))
addFactor!(fg, [:x2;:x3], pp, graphinit=false)

##

# nlsolve failing ?
meas = [10.0, 0.0, 1.0471975511965976]
X1 = [15.000000000016204, 8.660254037814505, 2.0943951023931953]
# X2 = [0.0004014798555661571, 0.0006953833518543132, 1.206557444798529e-10]
X2 = [10.00004891350537; 17.320479835550103; 4.498439149584132e-6]

M = getManifold(Pose2())
ϵ = getPointIdentity(Pose2())

X = hat(M, ϵ, meas) #measurement
p = exp(M, ϵ, hat(M, ϵ, X1)) # variable from
q = exp(M, ϵ, hat(M, ϵ, X2)) # variable to

## NLsolve not moving the cost value and "terminating too early"
# res = [3.469264352442895e-5, -2.002964655985515e-5, -3.1415887889926792]
# x = [9.999965307450084, 17.32052810543953, -3.864597114787283e-6]
# res = [3.469264352798169e-5, -2.0029646563407863e-5, -3.141588788992679]
# x = [9.999965307450083, 17.32052810543953, -3.864597114889992e-6]
# res = [3.469264352798169e-5, -2.0029646563407863e-5, -3.141588788992679]
# x = [9.999965307450083, 17.32052810543953, -3.864597114941346e-6]

# @test isapprox( 

resVal = calcFactorResidualTemporary( pp, (Pose2, Pose2), X, (p, q))

@test isapprox( resVal[1:2], [0;0], atol=1e-4)
@test isapprox( abs(resVal[3]), π, atol=1e-4)

# test positive angle
X2 = [10.00004891350537; 17.320479835550103; π]
q = exp(M, ϵ, hat(M, ϵ, X2)) # variable to
resVal = calcFactorResidualTemporary( pp, (Pose2, Pose2), X, (p, q))
@test isapprox( resVal[1:2], [0;0], atol=1e-4)
@test isapprox( abs(resVal[3]), 0, atol=1e-4)

# test negative angle
X2 = [10.00004891350537; 17.320479835550103; -π]
q = exp(M, ϵ, hat(M, ϵ, X2)) # variable to
resVal = calcFactorResidualTemporary( pp, (Pose2, Pose2), X, (p, q))
@test isapprox( resVal[1:2], [0;0], atol=1e-4)
@test isapprox( abs(resVal[3]), 0, atol=1e-4)


##

refVal = accumulateFactorMeans(fg, [:x2f1; :x2x3f1])

@test isapprox(refVal[1:2], [10;17.32], atol=1e-2)
@test isapprox( abs(refVal[3]), pi, atol = 1e-2)


##

end


@testset "Test canonical FG Honeycomb generation" begin

##

# build the graph
fg = RoME.generateCanonicalFG_Honeycomb!()

# check that pose :x3 has rotation near +-pi
t,m,gn = IIF._checkVariableByReference(fg, :x2, r"x\d+", Pose2, getFactorType(fg, :x2x3f1),destPrefix=:x, srcNumber=2)

@test isapprox( m.suggested[1], 10, atol = 1e-1)
@test isapprox( m.suggested[2], 17.32, atol = 1e-1)

@test isapprox( abs(m.suggested[3]), π, atol = 1e-2)


##

end

#