# basic Pose2 convolution test

using RoME
using Test

##

@testset "test basic Pose2 convolution based on MvNormal" begin
##

fg = initfg()

addVariable!(fg, :x0, Pose2)
u0 = ArrayPartition([0;0.0], [1 0; 0 1.])
M = getManifold(Pose2) # TODO add better dispatch to simplify
x0 = [AMP.makePointFromCoords(M,0.01*randn(3),u0) for _ in 1:100];
X0 = manikde!(M, x0)
initVariable!(fg, :x0, X0)

addVariable!(fg, :x1, Pose2)
addFactor!(fg, [:x0;:x1], Pose2Pose2(MvNormal([10;0;pi], 0.1*diagm([1;1;1]))), graphinit=false)

## now do a convolution

X1_ = approxConvBelief(fg, :x0x1f1, :x1)
x1_ = getPoints(X1_) .|> x->AMP.makeCoordsFromPoint(M,x)

μ_x1_ = mean(X1_)
Σ_x1_ =  cov(M, getPoints(X1_))

@test isapprox( submanifold_component(μ_x1_,1), [10; 0], atol=0.2)
@test isapprox( submanifold_component(μ_x1_,2), [-1 0; 0 -1], atol=0.2)

@test isapprox( Σ_x1_, [0.15 0 0; 0 0.15 0; 0 0 0.15], atol=0.4)

# check X1_ rotation's coordinates is spread between +pi and -pi

ψ_ = (x->x[3]).(x1_)

# bi-modality should occur since ±pi is the same value on-manifold
@test 20 < sum(2 .< ψ_)
@test 20 < sum(ψ_ .< -2)
# numerical check
@test sum(-2 .< ψ_ .< 2) < 1
# on-manifold check
@test sum(3.15 .< abs.(ψ_)) < 1


## now test the deconvolution

pts, meas = approxDeconv(fg, :x0x1f1)
X12_ = manikde!(M, pts)

# check that deconv is good

@test mmd(M, pts, meas) < 0.001

##
end





