# use Pose2 and test that the manifold products are working
using Statistics
using Test
using RoME, Distributions
# using KernelDensityEstimate

# using Gadfly
# using RoMEPlotting

using ApproxManifoldProducts
# const AMP = ApproxManifoldProducts

##


@testset "Testing hand constructed product of SE(2) equivalent..." begin

## Reuse variables

pt = Pose2()

##



# Pose densities on SE3 equivalent
pts1a = 0.1*randn(1,100)
pts1b = 0.1*randn(1,100)
pts1c = TransformUtils.wrapRad.( 0.1*randn(1,100) )
p = manikde!([pts1a;pts1b;pts1c], Pose2().manifolds )


#

pts2a = 3.0*randn(1,100).+4.0
pts2b = 3.0*randn(1,100)
pts2c = TransformUtils.wrapRad.(0.15*randn(1,100).+0.7pi)
q = manikde!([pts2a;pts2b;pts2c], Pose2().manifolds  )




pts3a = 3.0*randn(1,100).-4.0
pts3b = 3.0*randn(1,100)
pts3c = TransformUtils.wrapRad.(0.15*randn(1,100).-0.7pi)
r = manikde!([pts3a;pts3b;pts3c], Pose2().manifolds  )




# pl = plotPose(pt, [q; r], levels=1)

qr = manifoldProduct([q;r], pt.manifolds, Niter=3)

# pl = plotPose(pt, [q; r; qr], levels=2, c=["cyan";"cyan";"red"])


##




# Pose densities on SE3 equivalent
# pts1a = 0.1*randn(1,100)
# pts1b = 0.1*randn(1,100)
# pts1c = TransformUtils.wrapRad.( 0.1*randn(1,100) )
# p = manikde!([pts1a;pts1b;pts1c], Pose2().manifolds )


R = [1 1; -1 1]./sqrt(2)
I2 = [16.0 0; 0 0.5]
mvn = MvNormal([0.0;0.0],  R*I2*(R'))
pts1 = rand(mvn, 100)
# Gadfly.plot(x=pts1[1,:],y=pts1[2,:],Geom.hexbin )

pts2a = pts1[1:1,:].+5.0
pts2b = pts1[2:2,:].+5.0
pts2c = TransformUtils.wrapRad.(0.25*randn(1,100).+pi/2)
q = manikde!([pts2a;pts2b;pts2c], Pose2().manifolds  )


R = [1 -1; 1 1]./sqrt(2)
I2 = [16.0 0; 0 0.5]
mvn = MvNormal([0.0;0.0],  R*I2*(R'))
pts2 = rand(mvn, 100)
# Gadfly.plot(x=pts2[1,:],y=pts2[2,:],Geom.hexbin )

pts3a = pts2[1:1,:]
pts3b = pts2[2:2,:]
pts3c = TransformUtils.wrapRad.(0.25*randn(1,100).-pi/2)
r = manikde!([pts3a;pts3b;pts3c], Pose2().manifolds  )


##

qr = manifoldProduct([q;r], pt.manifolds, Niter=8)

# pl = plotPose(pt, [q; r; qr], levels=3, c=["cyan";"cyan";"red"])


##


# Pose densities on SE3 equivalent
pts1a = 0.5*randn(1,100)
pts1b = 0.5*randn(1,100)
pts1c = TransformUtils.wrapRad.( 0.1*randn(1,100) )
p = manikde!([pts1a;pts1b;pts1c], Pose2().manifolds )


pts2a = 0.5*randn(1,100)
pts2b = 1.0*randn(1,100)
pts2c = TransformUtils.wrapRad.(0.15*randn(1,100).+2pi/3)
q = manikde!([pts2a;pts2b;pts2c], Pose2().manifolds  )


pts3a = 1.0*randn(1,100)
pts3b = 0.5*randn(1,100)
pts3c = TransformUtils.wrapRad.(0.15*randn(1,100).-2pi/3)
r = manikde!([pts3a;pts3b;pts3c], Pose2().manifolds  )



pqr = manifoldProduct([p;q;r], pt.manifolds)

# pl = plotPose(pt, [p; q; r; pqr], levels=1, c=["cyan";"cyan";"cyan";"red"])


end # testset



## Save the latest plot to file


# pl |> SVG("/tmp/test.svg",12cm,10cm)
# @async run(`eog /tmp/test.svg`)
#
#
# pl |> PNG("/tmp/test.png",12cm,10cm)
# @async run(`eog /tmp/test.png`)


## repeat one of the examples but in a factor graph


@testset "Factor graph based test of SE(2) equivalent using only Priors..." begin

N = 100

#test as priors on one node

pts2a = 3.0*randn(1,N)
pts2b = 3.0*randn(1,N)
pts2c = TransformUtils.wrapRad.(0.15*randn(1,N).+0.9pi)
q = manikde!([pts2a;pts2b;pts2c], Pose2().manifolds  )


pts3a = 3.0*randn(1,N).+4.0
pts3b = 3.0*randn(1,N).+4.0
pts3c = TransformUtils.wrapRad.(0.15*randn(1,N).-0.9pi)
r = manikde!([pts3a;pts3b;pts3c], Pose2().manifolds  )


fg = initfg()

addVariable!(fg, :x1, Pose2)

addFactor!(fg, [:x1], Prior(q))
addFactor!(fg, [:x1], Prior(r))

# writeGraphPdf(fg, show=true)


stuff = IIF.localProduct(fg, :x1)

@test 0.90*N < sum( -10 .< getPoints(stuff[1])[1,:] .< 15)
@test 0.90*N < sum( -10 .< getPoints(stuff[1])[2,:] .< 15)

@test 0.90*N < sum(2.7 .< abs.(getPoints(stuff[1])[3,:]) .< 1.01*pi)


##
# plotPose(Pose2(), [getVertKDE(fg, :x1)])
#
# plotPose(Pose2(), [stuff[2];stuff[1]], levels=2, c=["cyan";"cyan";"red"])

##

end







@testset "Factor graph based test of SE(2) equivalent using Priors and Conditional Likelihoods..." begin

N = 100

@warn "not implemented yet"


end



##
