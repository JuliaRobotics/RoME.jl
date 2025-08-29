using KernelDensityEstimate, RoME, Gadfly



XY = Float64[]
for i in 1:200  XY=[XY;pol2cart([100.0+randn();2pi*rand()-pi],ones(2))[1]]; end
xy1 = reshape(XY[1:200],2,100);
xy2 = reshape(XY[201:400],2,100);
xy1[1,:] -= 25.0
xy2[1,:] += 25.0


plot(layer(x=xy1[1,:],y=xy1[2,:],Geom.point),layer(x=xy2[1,:],y=xy2[2,:],Geom.point))


p1 = kde!(xy1);
p2 = kde!(xy2);
p12 = p1*p2;

plotKDE(marginal(p12,[1]))
plotKDE(marginal(p12,[2]))
