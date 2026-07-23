

using LieGroups, Manifolds
using RoME
import Rotations as Rot_
using StaticArrays
using LinearAlgebra
using IncrementalInference.Manopt

## variable options with different tangent spaces are different
#
# 1. LeftInvariantMetricSE(2)
# 2. SpecialEuclideanGroup(2; variant = :right)
# 3. TranslationGroup(2) × SpecialOrthogonalGroup(2)

p = ArrayPartition(SA[1, 2.], SMatrix{2, 2}(Rot_.RotMatrix2(0.1)))
q = ArrayPartition(SA[3, 5.], SMatrix{2, 2}(Rot_.RotMatrix2(0.2)))

M1 = RoME.LeftInvariantMetricSE(2) # why SE in name here? 
M2 = SpecialEuclideanGroup(2; variant = :right)
M3 = TranslationGroup(2) × SpecialOrthogonalGroup(2)

X1 = log(M1, p, q)
X2 = log(M2, p, q)
X3 = log(M3, p, q)


DFG.@defStateType(
    Pose_SE2,
    SpecialEuclideanGroup(2; variant = :right),
    ArrayPartition(SA[0; 0.0], SA[1 0; 0 1.0])
)

DFG.@defStateType(
    Pose_Tr2xSO2,
    TranslationGroup(2) × SpecialOrthogonalGroup(2),
    ArrayPartition(SA[0; 0.0], SA[1 0; 0 1.0])
)


DFG.@defObservationType(
    Prior_SE2,
    PriorObservation,
    SpecialEuclideanGroup(2; variant = :right)
)
DFG.@defObservationType(
    Prior_Tr2xSO2,
    PriorObservation,
    TranslationGroup(2) × SpecialOrthogonalGroup(2)
)

function (cf::CalcFactor{<:Prior_SE2})(m, p)
    M = getManifold(Prior_SE2)
    X = log(M, p, m)
    return vee(LieAlgebra(M), X)
end
function (cf::CalcFactor{<:Prior_Tr2xSO2})(m, p)
    M = getManifold(Prior_Tr2xSO2)
    X = log(M, p, m)
    return vee(LieAlgebra(M), X)
end


## ========================= KINEMATIC MECHANISM  --  i.e. non-smooth path through variable epoch
## hex pose velocity points along hex



DFG.@defObservationType(
    TestKinP2P2,
    RelativeObservation,
    RoME.LeftInvariantMetricSE(2)
)
function (cf::CalcFactor{<:TestKinP2P2})(X, p, q)
    M = getManifold(TestKinP2P2)
    X̂ = log(M, p, q)
    return vee(LieAlgebra(M), X - X̂)
end





fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :w_P_b0, Pose_Tr2xSO2)
P = zeros(3)  # make on-manifold
prior = Prior_Tr2xSO2(MvNormal(P, 0.01 * Matrix{Float64}(LinearAlgebra.I, 3, 3)))
addFactor!(fg, [:w_P_b0,], prior)

addVariable!(fg, :w_P_b1, Pose_Tr2xSO2)
addVariable!(fg, :w_P_b2, Pose_Tr2xSO2)
addVariable!(fg, :w_P_b3, Pose_Tr2xSO2)
addVariable!(fg, :w_P_b4, Pose_Tr2xSO2)
addVariable!(fg, :w_P_b5, Pose_Tr2xSO2)
addVariable!(fg, :w_P_b6, Pose_Tr2xSO2)

# addVariable!(fg, :l0, Position{2})

addFactor!(fg, [:w_P_b0, :w_P_b1], TestKinP2P2(MvNormal([10.0,0,pi/3], diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b1, :w_P_b2], TestKinP2P2(MvNormal([10.0,0,pi/3], diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b2, :w_P_b3], TestKinP2P2(MvNormal([10.0,0,pi/3], diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b3, :w_P_b4], TestKinP2P2(MvNormal([10.0,0,pi/3], diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b4, :w_P_b5], TestKinP2P2(MvNormal([10.0,0,pi/3], diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b5, :w_P_b6], TestKinP2P2(MvNormal([10.0,0,pi/3], diagm([0.1, 0.2, 0.01].^2)))) 


stopping_criterion=StopAfterIteration(1000) | StopWhenGradientNormLess(1e-4) | StopWhenStepsizeLess(1e-4)
debug = [:Iteration, :Cost," | ", :tension, " | ", :Stepsize, " | ", :GradientNorm, " | ", :last_step_successful, "\n", :Stop]

IIF.solveGraphParametric!(fg; init=false, stopping_criterion, debug)

mean(getBelief(fg, :w_P_b6, :parametric))





## ========================= DYNAMIC MECHANISM  --  i.e. smooth path through variable epoch
## hex pose velocity points tangent to circle path


DFG.@defObservationType TestDynP2P2 RelativeObservation LieGroups.SpecialEuclideanGroup(2; variant = :right)
function (cf::CalcFactor{<:TestDynP2P2})(X, p, q)
    M = getManifold(TestDynP2P2)
    X̂ = log(M, p, q)
    return vee(LieAlgebra(M), X - X̂)
end

p = RoME.getPointIdentity(SpecialEuclideanGroup(2; variant = :right))
q = ArrayPartition(SA[10, 0.], SMatrix{2, 2}(Rot_.RotMatrix2(pi/3)))

M = SpecialEuclideanGroup(2; variant = :right)
X = log(M, p, q)
Xc = vee(LieAlgebra(M), X)

fg = initfg()
getSolverParams(fg).graphinit = false

addVariable!(fg, :w_P_b0, Pose_SE2)
P = zeros(3)  # make on-manifold
prior = Prior_SE2(MvNormal(P, 0.01 * Matrix{Float64}(LinearAlgebra.I, 3, 3)))
addFactor!(fg, [:w_P_b0,], prior)


addVariable!(fg, :w_P_b1, Pose_SE2)
addVariable!(fg, :w_P_b2, Pose_SE2)
addVariable!(fg, :w_P_b3, Pose_SE2)
addVariable!(fg, :w_P_b4, Pose_SE2)
addVariable!(fg, :w_P_b5, Pose_SE2)
addVariable!(fg, :w_P_b6, Pose_SE2)

# addVariable!(fg, :l0, Position{2})

addFactor!(fg, [:w_P_b0, :w_P_b1], TestDynP2P2(MvNormal(Xc, diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b1, :w_P_b2], TestDynP2P2(MvNormal(Xc, diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b2, :w_P_b3], TestDynP2P2(MvNormal(Xc, diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b3, :w_P_b4], TestDynP2P2(MvNormal(Xc, diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b4, :w_P_b5], TestDynP2P2(MvNormal(Xc, diagm([0.1, 0.2, 0.01].^2)))) 
addFactor!(fg, [:w_P_b5, :w_P_b6], TestDynP2P2(MvNormal(Xc, diagm([0.1, 0.2, 0.01].^2)))) 


stopping_criterion=StopAfterIteration(1000) | StopWhenGradientNormLess(1e-4) | StopWhenStepsizeLess(1e-4)
debug = [:Iteration, :Cost," | ", :tension, " | ", :Stepsize, " | ", :GradientNorm, " | ", :last_step_successful, "\n", :Stop]

IIF.solveGraphParametric!(fg; init=false, stopping_criterion, debug)

mean(getBelief(fg, :w_P_b6, :parametric))


## end of main code, plotting section follows

if false
## plotting
    using GLMakie

    function plot2D(subfg, x_labels = sortDFG(listVariables(subfg; whereLabel = startswith("w_P"))))
        pnts = map(x_labels) do v
            val = mean(getBelief(fg, v, :parametric))
            Point2f(val.x[1][1:2])
        end
        fig = lines(pnts; axis = (aspect = DataAspect(),))
        θs = map(x_labels) do v
            R = mean(getBelief(subfg, v, :parametric)).x[2]
            atan(R[2, 1], R[1, 1])
        end
        
        scatter!(pnts; rotation = θs, markersize = 15, marker = '➤')

        return fig
    end

## plot the 2D trajectory of the poses

    fig = plot2D(fg)

    theta = range(0, 2π, length=100)
    circle_x = 5 .+ 10 .* cos.(theta)
    circle_y = 8.66 .+ 10 .* sin.(theta)
    lines!(circle_x, circle_y; color=:green)
    scatter!([5], [8.66]; color=:black, markersize=10)

    # lines!([0, 9.06] ./ 2, [0, -5.24] ./ 2; color=:red)

    fig
##
end

#