

# testing
using RoME, KernelDensityEstimate

meas = LinearRangeBearingElevation((3.0,3e-4),(0.2,3e-4))
# LinearRangeBearingElevation(Distributions.Normal{Float64}(μ=0.2, σ=0.001),Distributions.Normal{Float64}(μ=3.0, σ=0.001),Distributions.Uniform{Float64}(a=-0.25133, b=0.25133))
X, pts = 0.01*randn(6,200), zeros(3,200);
# X[1:3,:] = 0.1*randn(3,200)
for i in 1:200
	project!(meas, X, pts, i)
end
p1 = kde!(pts);

# X = rand(6,200)
for i in 1:200
	backprojectRandomized!(meas, pts, X, i)
end
p2 = kde!(X);




using Gadfly

axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
      plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )


axis = [[-2.5;2.5]';[-2.5;2.5]';[-2.5;2.5]';[-2pi;2pi]';[-2pi;2pi]';[-2pi;2pi]']
draw( PDF("/home/dehann/Desktop/test.pdf",30cm,20cm),
      plotKDE( p2, dimLbls=["x";"y";"z"; "roll"; "pitch"; "yaw"], axis=axis)  )
