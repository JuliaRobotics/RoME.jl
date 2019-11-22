# diff existing odo to delta-odo

# using Revise

using TransformUtils
using RoME
# using Gadfly
# Gadfly.set_default_plot_size(35cm, 25cm)

using Test


@testset "test odometry accumulation MutablePose2Pose2Gaussian..." begin

## load the existing odometry data
data = randn(3,1000);
data[1:2,:] .*= 0.1
data[3,:] .*= 0.01
cumdata = cumsum(data, dims=2)
dcs = deepcopy(cumdata[1:2,:])
cumdata[3,:] .= wrapRad.(cumdata[3,:])
cumdata[1:2,:] .= cumsum(cumdata[1:2,:], dims=2)

# Gadfly.plot(x=cumdata[1,:],y=cumdata[2,:], Geom.path)
# Gadfly.plot(y=cumdata[3,:], Geom.path)

# theta
ddcd = zeros(2,size(cumdata,2))
ddcd[1:2,2:end] = diff(cumdata[1:2,:], dims=2)
ddcd[1:2,1] = ddcd[1:2,2]
# Gadfly.plot(y=ddcd[1,:], Geom.path)
# Gadfly.plot(y=ddcd[2,:], Geom.path)
course = atan.(ddcd[2,:],ddcd[1,:])
course .= wrapRad.(course)
yaw = course# + cumdata[3,:]
yaw .= wrapRad.(yaw)
# Gadfly.plot(y=course, Geom.path)

XX = cumdata[1,:]
YY = cumdata[2,:]
TH = yaw


# Gadfly.plot(x=XX,y=YY, Geom.path)
# Gadfly.plot(y=TH, Geom.path)


## Draw the deltas to reintegrate later


dt = 1.0
DX = zeros(3,length(XX))
nXYT__ = zeros(3,size(DX,2))
nXYT__[:,1] = [XX[1];YY[1];TH[1]]
for i in 2:length(XX)
  wTbk = SE2([XX[i-1];YY[i-1];TH[i-1]])
  wTbk1 = SE2([XX[i];YY[i];TH[i]])
  bkTbk1 = wTbk\wTbk1
  DX[:,i] = se2vee(bkTbk1)

  # test
  nXYT__[:,i] .= se2vee(SE2(nXYT__[:,i-1])*SE2(DX[:,i]))
  # nXYT__[:,i] .= se2vee(SE2(nXYT__[:,i-1])*bkTbk1)
end

# Gadfly.plot(y=DX[1,:], Geom.path)
# Gadfly.plot(y=DX[2,:], Geom.path)
# Gadfly.plot(y=DX[3,:], Geom.path)


# Gadfly.plot(x=nXYT__[1,:],y=nXYT__[2,:], Geom.path)
# Gadfly.plot(y=nXYT__[3,:], Geom.path)


## Test accumulation function from OdometryUtils.jl


mpp = MutablePose2Pose2Gaussian(MvNormal([XX[1];YY[1];TH[1]], 1e-3*Matrix(LinearAlgebra.I, 3,3)))

nXYT = zeros(3,size(DX,2))
Qc = 1e-6*Matrix(LinearAlgebra.I, 3,3)
for i in 1:size(DX,2)
  RoME.accumulateDiscreteLocalFrame!(mpp,DX[:,i],Qc,dt)
  nXYT[:,i] .= mpp.Zij.μ
end


# Gadfly.plot(x=nXYT[1,:],y=nXYT[2,:], Geom.path)
# Gadfly.plot(y=nXYT[3,:], Geom.path)


@test norm(cumdata[1:2,:] - nXYT[1:2,:]) < 1e-6
@test norm(difftheta.(yaw, nXYT[3,:])) < 1e-3


end


#
