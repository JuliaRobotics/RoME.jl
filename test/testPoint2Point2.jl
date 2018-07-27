using RoME, Distributions


fg = initfg()

addNode!(fg, :x0, Point2)
addFactor!(fg, [:x0], PriorPoint2(MvNormal(zeros(2), eye(2))))

addNode!(fg, :x2, Point2)
addFactor!(fg, [:x0;:x2], Point2Point2(MvNormal([10;0.0], eye(2))))


tree = wipeBuildNewTree!(fg)
