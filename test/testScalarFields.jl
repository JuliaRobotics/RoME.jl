# test scalar field

# using Revise
using Test
using ImageCore, ImageIO
using TensorCast
using Interpolations
using RoME

import IncrementalInference: LevelSetGridNormal

##

@testset "Basic low-res ScalarField localization" begin
##

# # load dem (18x18km span, ~17m/px)
x_min, x_max = -9000, 9000
y_min, y_max = -9000, 9000
# north is regular map image up
global img
x, y, img = RoME.generateField_CanyonDEM(1, 100, x_is_north=false, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)



## modify to generate elevation measurements (data/smallData as in Boxy) and priors

dem = Interpolations.LinearInterpolation((x,y), img) # interpolated DEM
elevation(p) = dem(getPPE(fg, p, :simulated).suggested[1:2]'...)
sigma_e = 0.01 # elevation measurement uncertainty

## test buildDEMSimulated to ensure interpolation matches raw data 
im = (j->((i->dem(i,j)).(x))).(y);
@cast im_[i,j] := im[j][i];
@test norm(im_ - img) < 1e-10


##


function cb(fg_, lastpose)
  global dem, img
  
  # query DEM at ground truth
  z_e = elevation(lastpose)
  
  # generate noisy measurement
  @info "Callback for DEM LevelSet priors" lastpose ls(fg_, lastpose) z_e
  
  # create prior
  hmd = LevelSetGridNormal(img, (x,y), z_e, sigma_e, N=10000, sigma_scale=1)
  pr = PartialPriorPassThrough(hmd, (1,2))
  addFactor!(fg_, [lastpose], pr, tags=[:DEM;], graphinit=false, nullhypo=0.1)
  nothing
end


## Testing 

# 0. init empty FG w/ datastore
fg = initfg()
# ensure specific solve settings
getSolverParams(fg).useMsgLikelihoods = true
getSolverParams(fg).graphinit = false
# getSolverParams(fg).treeinit = true

storeDir = joinLogPath(fg,"data")
mkpath(storeDir)
datastore = FolderStore{Vector{UInt8}}(:default_folder_store, storeDir) 
addBlobStore!(fg, datastore)

# new feature, going to temporarily disable as WIP
getSolverParams(fg).attemptGradients = false

##

# 1. load DEM into the factor graph
# point uncertainty - 2.5m horizontal, 1m vertical
# horizontal uncertainty chosen so that 3sigma is approx half the resolution
if false
  sigma = diagm([2.5, 2.5, 1.0])
  @time loadDEM!(fg, img, (x), (y), meshEdgeSigma=sigma);
end

##

# 2. generate trajectory 

μ0 = [-7000;-2000.0;pi/2]
@time generateGraph_Helix2DSlew!(10, posesperturn=30, radius=1500, dfg=fg, μ0=μ0, postpose_cb=cb) # , graphinit=false , slew_x=1/20)
deleteFactor!(fg, :x0f1)



## optional prior at start

mu0 = getPPE(fg, :x0, :simulated).suggested
pr0 = PriorPose2(MvNormal(mu0, 0.01.*[1;1;1;]))
addFactor!(fg, [:x0], pr0)

##

tree = solveTree!(fg);

## check at least the first five poses

for lb in sortDFG(ls(fg,r"x\d+"))[1:4]
  sim = getPPE(fg, lb, :simulated).suggested
  ppe = getPPE(fg, lb).suggested
  @test isapprox(sim[1:2], ppe[1:2], atol=1000)
  @test isapprox(sim[3], ppe[3], atol=0.75)
end

##

@error "Skipping latter part of testScalarTest.jl, see #518"
# try

# for lb in sortDFG(ls(fg,r"x\d+"))[5:end]
#   sim = getPPE(fg, lb, :simulated).suggested
#   ppe = getPPE(fg, lb).suggested
#   @test isapprox(sim[1:2], ppe[1:2], atol=400)
#   @test isapprox(sim[3], ppe[3], atol=0.5)
# end

# catch
#   @error "ScalarField test failure on latter half poses"
# end

##
end

##

# using Cairo, RoMEPlotting
# Gadfly.set_default_plot_size(35cm,20cm)

# plotSLAM2D_KeyAndRef(fg)
# imshow(img)

##



