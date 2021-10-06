# ScalarField functions related to Images.jl

@info "Loading RoME.jl tools related to ImageCore.jl"

using .ImageCore

export generateField_CanyonDEM


"""
    $SIGNATURES

Loads a sample DEM (as if simulated) on a regular grid... It's the Grand Canyon, 18x18km, 17m
"""
function generateField_CanyonDEM( scale=1, N=100; 
                                  x_is_north=true,
                                  x_min::Real=-9000, x_max::Real=9000,
                                  y_min::Real=-9000, y_max::Real=9000)
  #
  filepath = joinpath(dirname(dirname(@__DIR__)), "data","CanyonDEM.png")
  img_ = load(filepath) .|> Gray
  img_ = scale.*Float64.(img_)
  
  N_ = minimum([N; size(img_)...])
  img = img_[1:N_, 1:N_]
  
  # flip image so x-axis in plot is North and y-axis is West (ie img[-north,west], because top left is 0,0)
  _img_ = if x_is_north
    _img = collect(img')
    reverse(_img, dims=2)
  else
    # flip so north is down along with Images.jl [i,j] --> (x,y)
    reverse(img_, dims=1)
  end
  
  x = range(x_min, x_max, length = size(_img_,1)) # North
  y = range(y_min, y_max, length = size(_img_,2)) # East
  
  return (x, y, _img_)
end



#