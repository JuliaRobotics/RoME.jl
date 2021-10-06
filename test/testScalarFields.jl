# test scalar field

using Test
using ImageCore
using RoME

##

@testset "Basic low-res ScalarField localization" begin
##

# # load dem (18x18km span, ~17m/px)
x_min, x_max = -9000, 9000
y_min, y_max = -9000, 9000
# north is regular map image up
x, y, img_ = RoME.generateField_CanyonDEM(1, 100, x_is_north=false, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
# flip so north is down along with Images.jl [i,j] --> (x,y)
img = reverse(img_, dims=1)







##
end





##