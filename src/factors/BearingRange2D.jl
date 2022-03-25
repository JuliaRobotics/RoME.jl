
#-------------------------------------------------------------------------------
# bearing and range available

"""
    $TYPEDEF

Bearing and Range constraint from a Pose2 to Point2 variable.
"""
mutable struct Pose2Point2BearingRange{B <: IIF.SamplableBelief, R <: IIF.SamplableBelief} <: IIF.AbstractManifoldMinimize
  bearing::B
  range::R
end

getManifold(::IIF.InstanceType{<:Pose2Point2BearingRange}) = ProductGroup(ProductManifold(SpecialOrthogonal(2), TranslationGroup(1)))

# What a hot mess, TODO, clean up
function preambleCache(dfg::AbstractDFG, vars::AbstractVector{<:DFGVariable}, fct::Pose2Point2BearingRange)
  Mz = getManifold(fct)
  e0 = ProductRepr([1 0; 0 1.], [0.])
  Mp = SpecialEuclidean(2)
  Mp_0 = identity_element(Mp)
  T0_c = vee(Mp, Mp_0, Mp_0)

  smpl = hat(Mz, e0, [rand(fct.bearing), rand(fct.range)])
  bZ_c = vee(Mz, e0, smpl) # Manifolds.Identity(Mz)
  return (;Mz, e0, Mp, Mp_0, T0_c, bZ_c)
end

function getSample(cf::CalcFactor{<:Pose2Point2BearingRange})
  # defaults, TODO better reuse
  M = cf.cache.Mz
  e0 = cf.cache.e0

  # vector of tangents
  smpl = hat(M, e0, [rand(cf.factor.bearing), rand(cf.factor.range)])

  # return IIF `::Tuple` format
  return smpl
end


function IIF.getMeasurementParametric(s::Pose2Point2BearingRange{<:Normal, <:Normal})

  meas = [mean(s.bearing), mean(s.range)]
  iΣ = [1/var(s.bearing)             0;
                      0  1/var(s.range)]

  return meas, iΣ
end

function (cfo::CalcFactor{<:Pose2Point2BearingRange})(meas, wP, wL)
  Mp = cfo.cache.Mp
  Mt = Mp.manifold.manifolds[1]
  Mr = Mp.manifold.manifolds[2]
  # excess
  Mp_0 = cfo.cache.Mp_0
  T0_c = cfo.cache.T0_c
  Mz = cfo.cache.Mz
  bZ_c = cfo.cache.bZ_c
  meas_θ = bZ_c[1]
  meas_r = bZ_c[2]
  
  
  #FIXME fix this factor to work correctly on manifolds
  δ_l, δθ = if false
    
    # compose rotations
    rRi = wP.parts[2]
    rRo = retract(Mr, rRi, meas.parts[1])
    
    # new Mp objects containing the rotation
    rTo = ProductRepr(wP.parts[1], rRo)
    oTl = ProductRepr([meas.parts[2];0], identity_element(Mr, rRi))
    
    # prediction landmark
    rTo = Manifolds.compose(Mp, rTo, oTl)
    
    # cartesian difference in predicted and estimated landmark
    δ_l = Manifolds.compose(Mt, inv(Mt, rTo.parts[1]), wL)
    
    # find the residuals as though in polar coordinates
    δθ = atan(δ_l[2], δ_l[1])
    
    δ_l, δθ
  else # FIXME: this is close to the old BR that worked, use to ensure tests are working
    
    vee!(Mp, T0_c, wP, log(Mp, Mp_0, wP))
    vee!(Mz, bZ_c, cfo.cache.e0, meas)
    # 1-bearing
    # 2-range
    # world frame
    
    θ = meas_θ + T0_c[3]
    mx = meas_r*cos(θ)
    my = meas_r*sin(θ)
    
    ex = wL[1] - (mx + T0_c[1])
    ey = wL[2] - (my + T0_c[2])
    
    δ_l = [ex; ey]
    
    δθ = atan((my + T0_c[2]), (mx + T0_c[1])) - atan(wL[2], wL[1])
    
    δ_l, δθ
  end
  
  δr = norm(δ_l)
  
  return [δθ, δr]
end
# quick check
# pose = (0,0,0),  bear = 0.0,  range = 10.0   ==>  lm = (10,0)
# pose = (0,0,0),  bear = pi/2,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = 0.0,  range = 10.0   ==>  lm = (0,10)
# pose = (0,0,pi/2),  bear = pi/2,  range = 10.0   ==>  lm = (-10,0)


# function (cfo::CalcFactor{<:Pose2Point2BearingRange})(measX, p, q)
#   #
#   M = SpecialEuclidean(2)
#   q_SE = ProductRepr(q, identity_element(SpecialOrthogonal(2), p.parts[2]))

#   X_se2 = log(M, identity_element(M, p), compose(M, inv(M, p), q_SE))
#   X = X_se2.parts[1]
#   # NOTE this: `X̂ = log(M, p, q_SE)` wrong for what we want
#   #go to tangent vector
#   return measX - X 
# end

# Support for database based solving

passTypeThrough(d::FunctionNodeData{Pose2Point2Range}) = d

Base.@kwdef struct PackedPose2Point2BearingRange <: AbstractPackedFactor
    bearstr::PackedSamplableBelief
    rangstr::PackedSamplableBelief
end

function convert( ::Type{PackedPose2Point2BearingRange}, d::Pose2Point2BearingRange )
  return PackedPose2Point2BearingRange( convert(PackedSamplableBelief, d.bearing), convert(PackedSamplableBelief, d.range) )
end

# TODO -- should not be resorting to string, consider specialized code for parametric distribution types and KDEs
function convert( ::Type{Pose2Point2BearingRange}, d::PackedPose2Point2BearingRange )
  Pose2Point2BearingRange( convert(SamplableBelief, d.bearstr), convert(SamplableBelief, d.rangstr) )
end
