using RoME
using Test


lm_prior_noise = 0.01
meas_noise = 0.25
odom_noise = 0.1
n_samples = 200

# initialize mean landmark locations
l0 = 0.0
l1 = 1.0
l2 = 4.0

# "Ground-truth" robot poses
x0 = 0.0
x1 = 1.0
x2 = 2.0
x3 = 4.0

# Initialize empty factor graph
fg = initfg()

# Place strong prior on locations of three "doors"
addVariable!(fg, Symbol("l0"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, Symbol("l1"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)))

addVariable!(fg, Symbol("l2"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l2], Prior(Normal(l2, lm_prior_noise)))

# Add first pose
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)

# Make first "door" measurement
addFactor!(fg, [:x0; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

# Add second pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

# # Gaussian transition model
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(1, odom_noise)))

# # Make second "door" measurement
addFactor!(fg, [:x1; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

# # Add third pose
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)

# # Gaussian transition model
addFactor!(fg, [:x1; :x2], LinearConditional(Normal(1, odom_noise)))

# Add fourth pose
# addVariable!(fg, :x3, ContinuousScalar, N=n_samples)

# Add odometry transition and new landmark sighting
# addFactor!(fg, [:x2, :x3], LinearConditional(Normal(2, odom_noise)))
# addFactor!(fg, [:x3; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])


tree = batchSolve!(fg)

# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)



## Plotting functions below

using RoMEPlotting

# ensureAllInitialized!(fg)
plotKDE(fg, [:x0;:x1;:x2])

spyCliqMat(tree, :l0)

spyCliqMat(tree, :x2)


#
