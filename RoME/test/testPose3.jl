# using Revise
using RoME
using LieGroups
using StaticArrays
using Test
using LinearAlgebra
using DistributedFactorGraphs

##
@testset "test SE(3) coordinates to homography and back" begin
    ##

    M = SpecialEuclideanGroup(3; variant = :right)

    C = 0.2 * randn(6)
    H = coordinates_to_homography(M, C)
    C_ = homography_to_coordinates(M, H)

    @test isapprox(C, C_)

    ##
end

##
@testset "Test Basic Pose3 :parametric and :default" begin
    fg = initfg()

    prior_distribution = PriorPose3(
        MvNormal([0.0, 0, 0, 0, -pi / 4, 0], diagm([0.1, 0.1, 0.1, 0.01, 0.01, 0.01]) .^ 2),
    )

    addVariable!(fg, :x0, Pose3)
    prior = addFactor!(fg, [:x0], prior_distribution)

    odo_distribution = MvNormal(
        [sqrt(2), 0, 0, 0, 0, pi / 2],
        diagm([0.1, 0.1, 0.1, 0.01, 0.01, 0.01] .^ 2),
    )

    for i = 1:4
        f = Symbol("x", i - 1)
        t = Symbol("x", i)
        addVariable!(fg, t, Pose3)
        addFactor!(fg, [f, t], Pose3Pose3(odo_distribution))
    end

    initAll!(fg)

    @time r = IIF.solveGraphParametric!(fg)

    M = getManifold(Pose3)

    p0 = DFG.refMeans(getState(fg, :x0, :parametric))[1]
    p1 = DFG.refMeans(getState(fg, :x1, :parametric))[1]
    p2 = DFG.refMeans(getState(fg, :x2, :parametric))[1]
    p3 = DFG.refMeans(getState(fg, :x3, :parametric))[1]
    p4 = DFG.refMeans(getState(fg, :x4, :parametric))[1]

    @test isapprox(M, p0, p4, atol = 0.001)

    np0 = mean(M, getVal(fg, :x0))
    np1 = mean(M, getVal(fg, :x1))
    np2 = mean(M, getVal(fg, :x2))
    np3 = mean(M, getVal(fg, :x3))
    np4 = mean(M, getVal(fg, :x4))

    @test isapprox(M, p0, np0, atol = 0.1)
    @test isapprox(M, p1, np1, atol = 0.1)
    @test isapprox(M, p2, np2, atol = 0.1)
    @test isapprox(M, p3, np3, atol = 0.1)
    @test isapprox(M, p4, np4, atol = 0.1)
end

@testset "Pose3Pose3RotOffset factor residual" begin
    # Direct numerical evaluation of the factor residual
    # Uses a pitch offset (rotation about y-axis) to create cross-coupling 
    # between translation and rotation components — if the adjoint transform
    # is wrong, numerical values will be incorrect.

    M = getManifold(RoME.Pose3Pose3RotOffset)
    SO3 = SpecialOrthogonalGroup(3)
    𝔰𝔬3 = LieAlgebra(SO3)
    𝔤 = LieAlgebra(M)

    # Rotation offset: sensor pitched -0.3 rad about body y-axis
    β = 0.3
    bRa = Matrix(exp(SO3, hat(𝔰𝔬3, SA[0.0, -β, 0.0])))
    aRb = transpose(bRa)  # Ry(+0.3)

    # x0 at identity, x1 displaced 2m forward + 0.2 rad yaw in body frame
    x0 = ArrayPartition(SA[0.0, 0, 0], SA[1.0 0 0; 0 1 0; 0 0 1])
    # For LeftInvariantMetricSE (product geodesics):
    # body-frame translation is R_p^T*(t_q - t_p) and rotation is log_SO3(R_p^T*R_q)
    # With x0 at identity: body_trans = t_x1, body_rot = log_SO3(R_x1)
    x1 = ArrayPartition(SA[2.0, 0, 0], Matrix(exp(SO3, hat(𝔰𝔬3, SA[0.0, 0, 0.2]))))

    # Body-frame relative motion (what log(M, x0, x1) should give)
    body_coords = SA[2.0, 0.0, 0.0, 0.0, 0.0, 0.2]
    bX = hat(𝔤, body_coords, ArrayPartition)

    # Compute expected sensor-frame measurement via adjoint
    # Ad_{(0, aRb)}(v, ω) = (aRb*v, aRb*ω) since translation is zero
    # aTb = ArrayPartition(SA[0.0, 0, 0], aRb)
    # aX = adjoint(M, aTb, bX)

    bTa = ArrayPartition(SA[0.0, 0, 0], bRa)
    aX = diff_left_compose(base_lie_group(M), x0, bTa, bX)
    meas_coords = vee(𝔤, aX)

    # Verify analytical values: Ry(β) rotates body vectors into sensor frame
    # Translation: Ry(0.3) * [2, 0, 0] = [2cos(β), 0, -2sin(β)]
    @test isapprox(meas_coords[1], 2 * cos(β), atol = 1e-10)
    @test isapprox(meas_coords[2], 0.0, atol = 1e-10)
    @test isapprox(meas_coords[3], -2 * sin(β), atol = 1e-10)
    # Rotation: Ry(0.3) * [0, 0, 0.2] = [0.2*sin(β), 0, 0.2*cos(β)]
    @test isapprox(meas_coords[4], 0.2 * sin(β), atol = 1e-10)
    @test isapprox(meas_coords[5], 0.0, atol = 1e-10)
    @test isapprox(meas_coords[6], 0.2 * cos(β), atol = 1e-10)

    # Factor should give zero residual with correct bRa
    obs = RoME.Pose3Pose3RotOffset(MvNormal(Vector(meas_coords), 0.01 * I(6)))
    res = calcFactorResidualTemporary(obs, (Pose3, Pose3, RoME.Rotation3), aX, (x0, x1, bRa))
    @test norm(res) < 1e-10

    # Wrong bRa (arbitrary rotation) should give non-zero residual
    wrong_bRa = Matrix(exp(SO3, hat(𝔰𝔬3, SA[0.1, 0.2, -0.3])))
    res_wrong = calcFactorResidualTemporary(obs, (Pose3, Pose3, RoME.Rotation3), aX, (x0, x1, wrong_bRa))
    @test norm(res_wrong) > 0.1

    # Flipping bRa ↔ aRb should give non-zero residual
    res_flipped = calcFactorResidualTemporary(obs, (Pose3, Pose3, RoME.Rotation3), aX, (x0, x1, Matrix(aRb)))
    @test norm(res_flipped) > 0.1

    # --- Case 2: non-trivial poses (not at identity) ---
    # x0 facing +y (Rz(π/2)), x1 one step forward in body frame
    R0 = SA[0.0 -1 0; 1 0 0; 0 0 1]  # Rz(π/2)
    x0b = ArrayPartition(SA[1.0, 2, 0], R0)
    # body-frame trans: R0^T*(t1-t0) should be [1.5, 0, 0.5]; body rot: [0, 0.1, 0]
    body_trans = SA[1.5, 0.0, 0.5]
    body_rot = SA[0.0, 0.1, 0.0]
    body_coords_b = vcat(body_trans, body_rot)
    t1 = SA[1.0, 2, 0] + R0 * body_trans
    R1 = R0 * Matrix(exp(SO3, hat(𝔰𝔬3, body_rot)))
    x1b = ArrayPartition(t1, Matrix(R1))

    bX_b = hat(𝔤, body_coords_b, ArrayPartition)
    # aX_b = adjoint_action(M, aTb, bX_b)  # same bRa as before
    aX_b = diff_left_compose(base_lie_group(M), x0b, bTa, bX_b)

    obs_b = RoME.Pose3Pose3RotOffset(MvNormal(Vector(vee(𝔤, aX_b)), 0.01 * I(6)))
    res_b = calcFactorResidualTemporary(obs_b, (Pose3, Pose3, RoME.Rotation3), aX_b, (x0b, x1b, bRa))
    @test norm(res_b) < 1e-10
end

@testset "Test Basic Pose3 with Rotation offset :parametric and :default" begin
    # Solver integration test for Pose3Pose3RotOffset
    # Scenario: robot moves 1m forward per step (body x-axis), with 0.15 rad yaw per step.
    # Sensor has a z-axis rotation offset of α = 0.2 rad.
    # The measurement in sensor frame has rotated translation direction.

    fg = initfg()
    fg.solverParams.graphinit = false

    SO3 = SpecialOrthogonalGroup(3)
    𝔰𝔬3 = LieAlgebra(SO3)

    # Ground truth rotation offset
    α = 0.2
    δ = 0.15  # yaw change per step (body-frame rotation)
    bRa_true = Matrix(exp(SO3, hat(𝔰𝔬3, [0.0, 0, -α])))

    # Measurement in sensor frame:
    # aRb = Rz(+α) acts on body motion [1, 0, 0, 0, 0, δ]:
    #   translation: Rz(α)*[1,0,0] = [cos(α), sin(α), 0]
    #   rotation:    Rz(α)*[0,0,δ] = [0, 0, δ]  (z-axis invariant under Rz)
    odo_mean = SA[cos(α), sin(α), 0.0, 0.0, 0.0, δ]
    odo_distribution = MvNormal(odo_mean, diagm(SA[0.1, 0.1, 0.1, 0.01, 0.01, 0.01]) .^ 2)

    # Poses: x0 faces +y (Rz(π/2)), each step adds δ yaw and moves 1m in body x
    # Product-manifold log: body_trans = R_p^T*(t_q - t_p), so t_q - t_p = R_p*[1,0,0]
    # x0: [0, 0, 0], Rz(π/2)
    # x1: [0, 0, 0] + Rz(π/2)*[1,0,0] = [0, 1, 0], Rz(π/2 + δ)
    # x2: [0, 1, 0] + Rz(π/2+δ)*[1,0,0] = [-sin(δ), 1+cos(δ), 0], Rz(π/2 + 2δ)

    addVariable!(fg, :x0, Pose3)
    addFactor!(
        fg,
        [:x0],
        PriorPose3(
            MvNormal(
                SA[0.0, 0, 0, 0, 0, pi / 2],
                diagm(SA[0.1, 0.1, 0.1, 0.01, 0.01, 0.01]) .^ 2,
            ),
        ),
    )

    addVariable!(fg, :bRa, RoME.Rotation3)

    for i = 1:2
        f = Symbol("x", i - 1)
        t = Symbol("x", i)
        addVariable!(fg, t, Pose3)
        addFactor!(fg, [f, t, :bRa], RoME.Pose3Pose3RotOffset(odo_distribution))
    end

    # Prior on x2 at analytically correct position
    x2_pos = SA[-sin(δ), 1 + cos(δ), 0.0]
    x2_yaw = pi / 2 + 2δ
    addFactor!(
        fg,
        [:x2],
        PriorPose3(
            MvNormal(
                SA[x2_pos[1], x2_pos[2], x2_pos[3], 0.0, 0.0, x2_yaw],
                diagm(SA[0.1, 0.1, 0.1, 0.01, 0.01, 0.01]) .^ 2,
            ),
        ),
    )

    #FIXME needs IIF autoinitParametric! update PR #TBD
    # IIF.autoinitParametric!(fg)
    r = IIF.solveGraphParametric!(fg; init = false)

    M = getManifold(Pose3)

    p0 = DFG.refMeans(getState(fg, :x0, :parametric))[1]
    p1 = DFG.refMeans(getState(fg, :x1, :parametric))[1]
    p2 = DFG.refMeans(getState(fg, :x2, :parametric))[1]

    # Expected poses
    R_x0 = SA[0.0 -1 0; 1 0 0; 0 0 1]  # Rz(π/2)
    R_x1 = Matrix(exp(SO3, hat(𝔰𝔬3, [0.0, 0, pi / 2 + δ])))
    R_x2 = Matrix(exp(SO3, hat(𝔰𝔬3, [0.0, 0, x2_yaw])))

    @test isapprox(M, p0, ArrayPartition([0, 0.0, 0], R_x0), atol = 1e-3)
    @test isapprox(M, p1, ArrayPartition([0, 1.0, 0], R_x1), atol = 1e-3)
    @test isapprox(M, p2, ArrayPartition(Vector(x2_pos), R_x2), atol = 1e-3)

    # bRa should recover the z-rotation offset
    bRa_est = DFG.refMeans(getState(fg, :bRa, :parametric))[1]
    @test isapprox(SO3, bRa_est, bRa_true, atol = 1e-3)

    # Non-parametric: bRa cannot be initialized automatically (no prior on it)
    initAll!(fg)
    @test_broken isInitialized(fg, :bRa)

    # Initialize bRa manually near the true value and solve
    initVariable!(fg, :bRa, MvNormal([0.0, 0, -α], diagm([0.1, 0.1, 0.1]) .^ 2))
    solveGraph!(fg)

    np0 = mean(getBelief(fg, :x0))
    np1 = mean(getBelief(fg, :x1))
    np2 = mean(getBelief(fg, :x2))
    @test isapprox(M, np0, ArrayPartition([0, 0.0, 0], R_x0), atol = 2e-1)
    @test isapprox(M, np1, ArrayPartition([0, 1.0, 0], R_x1), atol = 2e-1)
    @test isapprox(M, np2, ArrayPartition(Vector(x2_pos), R_x2), atol = 2e-1)
    @test isapprox(IIF.calcMeanMaxSuggested(fg, :bRa).suggested, [0, 0, -α], atol = 2e-1)
end
