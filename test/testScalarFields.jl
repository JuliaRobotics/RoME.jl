# test scalar field

# using Revise
using Test
using ImageCore, ImageIO, ImageFiltering
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

# remember for Github CI to use something near N=10000 (24Q1)
global N = 100000

# https://docs.juliadsp.org/stable/convolutions/

## modify to generate elevation measurements (data/smallData as in Boxy) and priors

# elevation model
dem = Interpolations.LinearInterpolation((x,y), img) # interpolated DEM
elevation(p) = dem(getPPE(fg, p, :simulated).suggested[1:2]'...)
sigma_e = 0.01 # elevation measurement uncertainty
@test 0.46 < dem(0.0, 0.0) < 0.48

## test buildDEMSimulated to ensure interpolation matches raw data 
im = (j->((i->dem(i,j)).(x))).(y);
@cast im_[i,j] := im[j][i];
@test norm(im_ - img) < 1e-10

# magnitude of slope model
# dz/dx = z(x_{k+1}) - z(x_{k-1}) / (x_{k+1} - x_{k-1})
dimg = imgradients(img, KernelFactors.ando3) # dimg = gx, gy
mag = sqrt.(dimg[1].^2 .+ dimg[2].^2)        # gradient magnitude 
dsm = Interpolations.LinearInterpolation((x,y), mag)    # interpolated DEM
slope(p) = dsm(getPPE(fg, p, :simulated).suggested[1:2]'...)
sigma_m = 1e-3  # slope measurement uncertainty (based on slopes under 6e-2)
@test 0 < dsm(0.0, 0.0) < 0.07

# direction of slope as cos(grad) and sin(grad)
gimg = complex.(dimg[1], dimg[2])
gimgh = gimg ./ norm.(gimg)
gimg_cos = real.(gimgh)
gimg_sin = imag.(gimgh)
dsc = Interpolations.LinearInterpolation((x,y), gimg_cos)    # interpolated cos(grad)
dss = Interpolations.LinearInterpolation((x,y), gimg_sin)    # interpolated sin(grad)
direction_cos(p) = dsc(getPPE(fg, p, :simulated).suggested[1:2]'...)
direction_sin(p) = dss(getPPE(fg, p, :simulated).suggested[1:2]'...)
sigma_s = 0.05


##

function postpose_cb(fg_, lastpose)
  global dem, img, N
  
  # query DEM at ground truth
  z_e = elevation(lastpose)
  z_m = slope(lastpose)
  z_cg = direction_cos(lastpose)
  z_sg = direction_sin(lastpose)
  
  # generate noisy measurement
  @info "Callback for DEM LevelSet priors" lastpose ls(fg_, lastpose) z_e
  
  # create prior
  lse = LevelSetGridNormal(img, (x,y), z_e, sigma_e; N, sigma_scale=1)
  lsm = LevelSetGridNormal(mag, (x,y), z_m, sigma_m; N, sigma_scale=1)  # elevation measurement
  ls_cg = LevelSetGridNormal(gimg_cos, (x,y), z_cg, sigma_s; N, sigma_scale=1)
  ls_sg = LevelSetGridNormal(gimg_sin, (x,y), z_sg, sigma_s; N, sigma_scale=1)

  hm_em_data = lse.heatmap.data .* lsm.heatmap.data .* ls_cg.heatmap.data .* ls_sg.heatmap.data
  σ_em = sqrt(sigma_e^2 + sigma_m^2 + sigma_s^2 + sigma_s^2)
  ls_em = LevelSetGridNormal(hm_em_data, (x,y), z_e*z_m*z_cg*z_sg, σ_em; N, sigma_scale=1)
  
  pr = PartialPriorPassThrough(ls_em, (1,2))
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
@time generateGraph_Helix2DSlew!(1; posesperturn=30, radius=1500, dfg=fg, postpose_cb, μ0) # , graphinit=false , slew_x=1/20)
# @time generateGraph_Helix2DSlew!(10; posesperturn=30, radius=1500, dfg=fg, postpose_cb, μ0) # , graphinit=false , slew_x=1/20)
# @time generateGraph_Helix2DSlew!(20; posesperturn=30, radius=1500, dfg=fg, postpose_cb, μ0) # , graphinit=false , slew_x=1/20)
deleteFactor!(fg, :x0f1)

## optional prior at start

mu0 = getPPE(fg, :x0, :simulated).suggested
pr0 = PriorPose2(MvNormal(mu0, 0.01.*[1;1;1;]))
addFactor!(fg, [:x0], pr0)


## EXPLORE MODEL

# just add, solve, repeat (set to max 6 for Github CI 24Q1)
for i in 41:60
  @info "Solving" i
  generateGraph_Helix2DSlew!(i; posesperturn=30, radius=1500, dfg=fg, postpose_cb, μ0) # , graphinit=false , slew_x=1/20)
  if 0 == (i % 2)
    solveGraph!(fg; multithread=true, storeOld=true);
  end
end

##

# tree = solveGraph!(fg; multithread=true, storeOld=true);
# [ solveGraph!(fg; multithread=true, storeOld=true) for _ in 1:5];

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

plotSLAM2D_KeyAndRef(fg)
# imshow(img)

##

# using ImageView

##






