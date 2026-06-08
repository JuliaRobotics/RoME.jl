using RoME
using LieGroups
using Test
using ForwardDiff

function propagate_with_sigmapoints(G_domain, G_codomain, p_mean, Xp_mean, Sigma_prior, Sigma_meas)
    N_p = size(Sigma_prior, 1)
    N_m = size(Sigma_meas, 1)
    N = N_p + N_m # Total dimensions = 6
    
    # Stack into a Joint Covariance Matrix
    Sigma_joint = zeros(N, N)
    Sigma_joint[1:N_p, 1:N_p] = Sigma_prior
    Sigma_joint[N_p+1:end, N_p+1:end] = Sigma_meas
    
    # Extract the principal axes of uncertainty via Cholesky
    L = cholesky(Sigma_joint).L
    
    # Standard Unscented Transform Weights (2N points)
    c = sqrt(N)
    weight = 1.0 / (2 * N)
    
    q_sigma_points = []
    
    # Generate and Propagate Sigma Points
    for i in 1:N
        for sign in [1.0, -1.0]
            # 6D perturbation vector
            delta = sign * c * L[:, i]
            
            delta_p = delta[1:N_p]
            delta_m = delta[N_p+1:end]
            
            # Apply perturbation to the prior pose p
            # (Mapping coordinates into the tangent space and taking the exponential)
            p_pert = exp(G_domain, p_mean, hat(G_domain, p_mean, delta_p))
            
            # Apply perturbation to the measurement
            X_meas_pert = Xp_mean + delta_m
            
            # Propagate to q using the factor's exact forward model
            # Since factor is: X_meas = log(M, p, q) => q = exp(M, p, X_meas)
            q_pert = exp(G_codomain, p_pert, hat(G_codomain, p_pert, X_meas_pert))
            
            push!(q_sigma_points, q_pert)
        end
    end
    
    # Reconstruct the Covariance at q
    # First, find the expected mean of q
    q_mean = exp(G_codomain, p_mean, hat(G_codomain, p_mean, Xp_mean))
    
    Sigma_q = zeros(N_p, N_p)
    for q_pert in q_sigma_points
        # Pull back each perturbed q into the tangent space of the mean
        diff_hat = log(G_domain, q_mean, q_pert)
        diff_coords = vee(G_domain, q_mean, diff_hat)
        
        # Standard sample covariance
        Sigma_q += weight * (diff_coords * diff_coords')
    end
    
    return Sigma_q
end

function propagate_with_jacobians(M_dom, M_cod, p, Xp_coords, Sigma_prior, Sigma_meas)
    
    q = exp(M_cod, p, hat(LieAlgebra(M_cod), Xp_coords, ArrayPartition))

    alg = LieAlgebra(M_cod)

    # residual function closure for Pose2Pose2 factor
    function residual(dp, dq, dX)
        p_pert = exp(M_dom, p, hat(alg, dp, ArrayPartition{eltype(dp)}))
        q_pert = exp(M_dom, q, hat(alg, dq, ArrayPartition{eltype(dq)}))
        Xp_pert = hat(alg, Xp_coords + dX, ArrayPartition{eltype(dX)})

        X̂p = log(M_cod, p_pert, q_pert)
        
        return vee(alg, Xp_pert - X̂p)
    end

    # $$\frac{\partial f}{\partial q} \frac{\partial q}{\partial X} + \frac{\partial f}{\partial X} = 0$$
    # Extract each Jacobian
    J_p = ForwardDiff.jacobian(dp -> residual(dp, zeros(3), zeros(3)), zeros(3))
    J_q = ForwardDiff.jacobian(dq -> residual(zeros(3), dq, zeros(3)), zeros(3))
    J_X = ForwardDiff.jacobian(dX -> residual(zeros(3), zeros(3), dX), zeros(3))

    q_J_p = -J_q \ J_p # How q changes when prior p changes
    q_J_X = -J_q \ J_X # How q changes when measurement X changes

    return q_J_p * Sigma_prior * q_J_p' + q_J_X * Sigma_meas * q_J_X'
end

DFG.@defObservationType ScrewPose2Pose2 RelativeObservation SpecialEuclideanGroup(2; variant = :right)
function (cf::CalcFactor{<:ScrewPose2Pose2})(X, p, q)
    M = getManifold(ScrewPose2Pose2)
    X̂ = log(M, p, q)
    return vee(M, p, X - X̂)
end


@testset "Propagate covariance on SE(2)" begin

    M = getManifold(ScrewPose2Pose2)
    p = getPointIdentity(M)
    
    # Xp_coords = [10.0, 0.0, pi/4] 
    Xp_coords = [10.0, 0.1, pi/8] 
    q = exp(M, p, hat(LieAlgebra(M), Xp_coords, ArrayPartition))
    
    # Define Covariances
    Sigma_prior = diagm([0.001, 0.002, 0.003])
    Sigma_meas  = diagm([0.03, 1.0, 0.01] .^ 2)

    # propagate
    M_dom = getManifold(Pose2)
    M_cod = getManifold(ScrewPose2Pose2)

    jac_cov_x2 = propagate_with_jacobians(M_dom, M_cod, p, Xp_coords, Sigma_prior, Sigma_meas)

    ut_cov_x2 = propagate_with_sigmapoints(M_dom, M_cod, p, Xp_coords, Sigma_prior, Sigma_meas)

    # --- Setup and Solve Graph ---
    fg = initfg()
    getSolverParams(fg).graphinit = false
    addVariable!(fg, :x1, Pose2)
    addVariable!(fg, :x2, Pose2)
    
    addFactor!(fg, [:x1], PriorPose2(MvNormal([0.0, 0.0, 0.0], Sigma_prior)))
    addFactor!(fg, [:x1; :x2], ScrewPose2Pose2(MvNormal(Xp_coords, Sigma_meas)))

    IIF.solveGraphParametric!(fg)
    
    x1 = getState(fg, :x1, :parametric)
    x2 = getState(fg, :x2, :parametric)

    # x1 should just match the prior
    @test isapprox(DFG.refMeans(x1)[1], p; atol = 1e-6)
    @test isapprox(DFG.refCovariances(x1)[1], Sigma_prior; atol = 1e-6)

    # x2 should match propagated mean and covariance
    @test isapprox(DFG.refCovariances(x2)[1], jac_cov_x2; atol = 1e-6)
    @test isapprox(DFG.refCovariances(x2)[1], ut_cov_x2; atol = 3e-3)

    @test isapprox(DFG.refMeans(x2)[1], q; atol = 1e-6)

end

@testset "Propagate covariance on LeftInvariantMetricSE(2)" begin

    M = getManifold(Pose2Pose2)
    p = getPointIdentity(M)
    
    # Xp_coords = [10.0, 0.0, pi/4] 
    Xp_coords = [10.0, 0.1, pi/8] 
    q = exp(M, p, hat(LieAlgebra(M), Xp_coords, ArrayPartition))
    
    # Define Covariances
    Sigma_prior = diagm([0.001, 0.002, 0.003])
    Sigma_meas  = diagm([0.03, 1.0, 0.01] .^ 2)

    # Propagate to x2
    M_dom = getManifold(Pose2)
    M_cod = getManifold(Pose2Pose2)

    jac_cov_x2 = propagate_with_jacobians(M_dom, M_cod, p, Xp_coords, Sigma_prior, Sigma_meas)

    ut_cov_x2 = propagate_with_sigmapoints(M_dom, M_cod, p, Xp_coords, Sigma_prior, Sigma_meas)

    # --- Setup and Solve Graph ---
    fg = initfg()
    getSolverParams(fg).graphinit = false
    addVariable!(fg, :x1, Pose2)
    addVariable!(fg, :x2, Pose2)
    
    addFactor!(fg, [:x1], PriorPose2(MvNormal([0.0, 0.0, 0.0], Sigma_prior)))
    addFactor!(fg, [:x1; :x2], Pose2Pose2(MvNormal(Xp_coords, Sigma_meas)))

    IIF.solveGraphParametric!(fg)
    
    x1 = getState(fg, :x1, :parametric)
    x2 = getState(fg, :x2, :parametric)

    # x1 should just match the prior
    @test isapprox(DFG.refMeans(x1)[1], p; atol = 1e-4)
    @test isapprox(DFG.refCovariances(x1)[1], Sigma_prior; atol = 1e-4)

    # x2 should match propagated mean and covariance
    @test isapprox(DFG.refCovariances(x2)[1], ut_cov_x2; atol = 3e-3)
    @test isapprox(DFG.refCovariances(x2)[1], jac_cov_x2; atol = 1e-4)
    @test isapprox(DFG.refMeans(x2)[1], q; atol = 1e-4)

end
