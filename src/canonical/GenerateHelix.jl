

## DEPRECATE BELOW TO GENERALIZED HELIX


"""
    $SIGNATURES

Generate an appoximate 2D helix (1 turn).

DevNotes
- TODO replace with generalized helix parameterized by a curve along z-axis  
"""
function _calcHelix2DApprox(; N_ppt::Integer = 20,
                              radius::Real = 0.5,
                              runback::Real = 5/7  )
  #
  @assert iseven(N_ppt) "N_ppt=$N_ppt must be an even number"

  # two part construction of helix, top and bottom
  theta_top = LinRange(0,pi,Int(N_ppt/2)) |> reverse
  theta_bot = - LinRange(0,pi,Int(N_ppt/2))
  # top half xy is easy, and get tangent gradient as angle from local
  top = exp.(im.*theta_top)
  angt = rem2pi.(angle.(top) .- pi/2, RoundNearest)
  # bottom half to be squashed in local x
  bot_ = exp.(im.*theta_bot)
  angb_ = rem2pi.(angle.(bot_) .- pi/2, RoundNearest)
  # and sqeeze on local x (not y)
  bx_ = (runback * (real.(bot_) .- 1) .+ 1)
  bot = bx_ .+ im.*imag.(bot_)
  # special care on squashed gradiens for bottom half 
  dydx = exp.(im.*angb_)
  dy = imag.(dydx)
  dx = runback .* real.(dydx)
  angb = angle.(dx .+ im.*dy)

  # scale to radius and offset for starting at local 0
  loop = radius .* (1 .+ vcat(top, bot[2:end]))
  ang = vcat(angt, angb[2:end])

  # return 2D array of data, rows are (x,y,theta) and columns are knot/pose points around helix2D
  return hcat(real.(loop), imag.(loop), ang)'
end


# TODO replace with generalized helix generator
function _calcHelix2DTurnsX(turns=1;
                            N_ppt::Integer = 20,  
                            radius::Real = 0.5,
                            runback::Real = 5/7  )
  #

  allpts = Matrix{Float64}(undef, 3, 0)
  for tn in 0:(turns-1)
    tmp_ = _calcHelix2DApprox(N_ppt=N_ppt, radius=radius, runback=runback)
    tmp_[1,:] .+= tn*(2radius*(1-runback))
    allpts = hcat(allpts, tmp_[:,2:end])
  end

  return allpts
end