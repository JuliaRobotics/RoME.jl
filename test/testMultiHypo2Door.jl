# load requried packages
using RoME
using Test


## parameters

lm_prior_noise = 0.01
meas_noise = 0.25
odom_noise = 0.1
n_samples = 100

# initialize mean landmark locations
l0 = 0.0
l1 = 10.0
l2 = 40.0

# "Ground-truth" robot poses
x0 = 0.0
x1 = 10.0
x2 = 20.0
x3 = 40.0

## Initialize empty factor graph
fg = initfg()

# Place strong prior on locations of three "doors"
addVariable!(fg, Symbol("l0"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l0], Prior(Normal(l0, lm_prior_noise)))

addVariable!(fg, Symbol("l1"), ContinuousScalar, N=n_samples)
addFactor!(fg, [:l1], Prior(Normal(l1, lm_prior_noise)))


# Add first pose
addVariable!(fg, :x0, ContinuousScalar, N=n_samples)

# Make first "door" measurement
addFactor!(fg, [:x0; :l0], LinearConditional(Normal(0, meas_noise)))
# addFactor!(fg, [:x0; :l0; :l1], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/2.0; 1.0/2.0])


# Add second pose
addVariable!(fg, :x1, ContinuousScalar, N=n_samples)

# Gaussian transition model
addFactor!(fg, [:x0; :x1], LinearConditional(Normal(x1-x0, odom_noise)))

# Make second "door" measurement
# addFactor!(fg, [:x1; :l1], LinearConditional(Normal(0, meas_noise)) )
addFactor!(fg, [:x1; :l0; :l1], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/2.0; 1.0/2.0])




## Add one more pose/odometry to invoke issue #236

# Add third pose
addVariable!(fg, :x2, ContinuousScalar, N=n_samples)
addFactor!(fg, [:x1; :x2], LinearConditional(Normal(x2-x1, odom_noise)))


# Add fourth pose
# addVariable!(fg, :x3, ContinuousScalar, N=n_samples)

# Add odometry transition and new landmark sighting
# addFactor!(fg, [:x2, :x3], LinearConditional(Normal(2, odom_noise)))
# addFactor!(fg, [:x3; :l0; :l1; :l2], LinearConditional(Normal(0, meas_noise)), multihypo=[1.0; 1.0/3.0; 1.0/3.0; 1.0/3.0])

## Do some debugging

writeGraphPdf(fg, show=true)

## Solve graph
tree = batchSolve!(fg)

# tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)



## Plotting functions below

using RoMEPlotting


ensureAllInitialized!(fg)
plotKDE(fg, [:x0;:x1])

plotKDE(fg, [:x0;:x1;:x2])

spyCliqMat(tree, :l0)

spyCliqMat(tree, :x2)



## specialized debugging


plotLocalProduct(fg, :x0)
plotLocalProduct(fg, :x1)


##

stuff = treeProductUp(fg, tree, :x2, :x2)


##

ensureAllInitialized!(fg)


tree = wipeBuildNewTree!(fg, drawpdf=true, show=true)
cliqorder = getCliqOrderUpSolve(tree)


spyCliqMat(cliqorder[end])



#
