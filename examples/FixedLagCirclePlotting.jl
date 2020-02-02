

function plotCirc10BA(fg_::AbstractDFG,
                      fg::AbstractDFG;
                      pointsList_::Vector{Symbol}=Symbol[],
                      pointsList::Vector{Symbol}=Symbol[],
                      filepath::AbstractString="/tmp/test.pdf",
                      lineColor_::String="gray40",
                      lineColor::String="green",
                      ellipseList_=ls(fg_),
                      ellipseColor_::String="gray60",
                      ellipseColor::String="gray20",
                      drawEllipse::Bool=true,
                      contourList_::Vector{Symbol}=Symbol[],
                      contourColor_::String="red",
                      contourLineWidth_=1pt,
                      levels::Int=3,
                      show::Bool=true,
                      width=10cm, height=8cm)
  #
  plfl4 = drawPosesLandms(fg_, spscale=1.5, manualColor=lineColor_, point_size=4pt, drawhist=false, contour=false, levels=2, lbls=false)

  for vsym in ellipseList_
    PL = plotCovEllipseLayer(fg_, vsym, points=vsym in pointsList_ || vsym in pointsList, ellipseColor=ellipseColor_)
    if vsym != :x10
      for ll in PL
        push!(plfl4.layers, ll)
      end
    end
  end

  if exists(fg, :x9)
    PL = plotCovEllipseLayer(fg, :x9, points=:x9 in pointsList, ellipseColor=ellipseColor, drawEllipse=drawEllipse)
    for ll in PL
      push!(plfl4.layers, ll)
    end
  end

  plfl5 = drawPosesLandms(fg, spscale=1.5, point_size=3pt, manualColor=lineColor, drawhist=false, contour=false, line_width=2pt)

  for ll in plfl4.layers
    push!(plfl5.layers, ll)
  end

  for conl in contourList_
    X = getKDE(fg_, conl)
    plcon = plotKDEContour(marginal(X,[1;2]), levels=levels, c=[contourColor_], line_width=contourLineWidth_)
    for ll in plcon.layers
      push!(plfl5.layers, ll)
    end
  end

  plfl5.coord = Coord.Cartesian(xmin=-20,xmax=30,ymin=-5,ymax=40)

  plfl5 |> PDF("$filepath", width, height);
  if show
    @async run(`evince $filepath`)
  end
  plfl5
end



#
