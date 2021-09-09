
"""
$(TYPEDEF)

Direction observation information of a `Point2` variable.
"""
mutable struct PriorPoint2{T <: IIF.SamplableBelief} <: IncrementalInference.AbstractPrior
  Z::T
  # empty constructor
  # PriorPoint2{T}() where T = new()
  # regular constructor
  # PriorPoint2{T}(dist::T) where {T <: IIF.SamplableBelief} = new{T}(dist)
end

# convenience helper and default object
PriorPoint2() = PriorPoint2(MvNormal(zeros(2),LinearAlgebra.diagm([0.01;0.01])))

DFG.getManifold(::PriorPoint2) = TranslationGroup(2)

function getSample(cfo::CalcFactor{<:PriorPoint2})
  return rand(cfo.factor.Z)
end

function (cf::CalcFactor{<:PriorPoint2})(meas, 	
                                        X1)	
  #	
  return meas[1:2] .- X1[1:2] 	
end

"""
$(TYPEDEF)
"""
mutable struct Point2Point2{D <: IIF.SamplableBelief} <: AbstractRelativeRoots
    Zij::D
    # empty constructor
    Point2Point2{T}() where T = new{T}()
    # regular constructor
    Point2Point2{T}(x::T) where {T <: IIF.SamplableBelief} = new{T}(x)
end
# convenience and default object helper
Point2Point2(x::T=MvNormal(zeros(2),LinearAlgebra.diagm([0.1;0.1]))) where {T <: IIF.SamplableBelief} = Point2Point2{T}(x)

function getSample(cfo::CalcFactor{<:Point2Point2})
  return rand(cfo.factor.Zij)
end
function (pp2r::CalcFactor{<:Point2Point2})(meas,
                                            xi,
                                            xj )
  #
  return meas[1:2] .- (xj[1:2] .- xi[1:2])
end




# ---------------------------------------------------------



"""
$(TYPEDEF)

Serialization type for `PriorPoint2`.
"""
mutable struct PackedPriorPoint2  <: IncrementalInference.PackedInferenceType
    str::String
    PackedPriorPoint2() = new()
    PackedPriorPoint2(x::String) = new(x)
end


function convert(::Type{PriorPoint2}, d::PackedPriorPoint2)
  # Cov = reshapeVec2Mat(d.vecCov, d.dimc)
  distr = convert(SamplableBelief, d.str)
  return PriorPoint2{typeof(distr)}(distr)
end
function convert(::Type{PackedPriorPoint2}, d::PriorPoint2)
  # v2 = d.mv.Σ.mat[:];
  return PackedPriorPoint2(convert(PackedSamplableBelief, d.Z))
end





"""
$(TYPEDEF)

Serialization type for `Point2Point2`.
"""
mutable struct PackedPoint2Point2 <: IncrementalInference.PackedInferenceType
    str::String
    PackedPoint2Point2() = new()
    PackedPoint2Point2(s::AS) where {AS <: AbstractString} = new(s)
end
function convert(::Type{Point2Point2}, d::PackedPoint2Point2)
  return Point2Point2( convert(SamplableBelief, d.str) )
end
function convert(::Type{PackedPoint2Point2}, d::Point2Point2)
  return PackedPoint2Point2( convert(PackedSamplableBelief, d.Zij) )
end

