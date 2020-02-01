

function plotCirc10BA(fg_::AbstractDFG,
                      fg::AbstractDFG;
                      pointsList::Vector{Symbol}=Symbol[],
                      filepath::AbstractString="/tmp/test.pdf",
                      show::Bool=true)
  #
  plfl4 = drawPosesLandms(fg_, spscale=1.5, manualColor="gray40", point_size=4pt, drawhist=false, contour=false, levels=2, lbls=false)

  for vsym in ls(fg)
    PL = plotCovEllipseLayer(fg_, vsym, points=vsym in pointsList, ellipseColor="gray60")
    if vsym != :x10
      for ll in PL
        push!(plfl4.layers, ll)
      end
    end
  end

  if exists(fg, :x9)
    PL = plotCovEllipseLayer(fg, :x9, points=false, ellipseColor="gray60")
    push!(plfl4.layers, PL[1])
  end

  plfl5 = drawPosesLandms(fg, spscale=1.5, point_size=3pt, manualColor="green", drawhist=false, contour=false)

  for ll in plfl4.layers
    push!(plfl5.layers, ll)
  end

  plfl5.coord = Coord.Cartesian(xmin=-20,xmax=30,ymin=-5,ymax=40)

  plfl5 |> PDF("$filepath", 10cm, 8cm);
  if show
    @async run(`evince $filepath`)
  end
  plfl5
end



#
