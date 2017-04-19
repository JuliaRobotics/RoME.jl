# inertial pose3

export
  PreintegralCompensationGradients,
  InertialPose3Container,
  oplus,
  ⊕,
  getSample,
  InertialPose3,
  PackedInertialPose3,
  PriorInertialPose3,
  PackedPriorInertialPose3,
  compare

abstract PreintContainer

type PreintegralCompensationGradients <: PreintContainer
  # First order terms
  dPdDa::Array{Float64,2}
  dVdDa::Array{Float64,2}
  dPdDw::Array{Float64,2}
  dVdDw::Array{Float64,2}
  dRdDw::Array{Float64,2}
  # second order
  # ...
  PreintegralCompensationGradients(;
        dPdDa::Array{Float64,2}=zeros(3,3),
        dVdDa::Array{Float64,2}=zeros(3,3),
        dPdDw::Array{Float64,2}=zeros(3,3),
        dVdDw::Array{Float64,2}=zeros(3,3),
        dRdDw::Array{Float64,2}=zeros(3,3)  ) = new(dPdDa,dVdDa,dPdDw,dVdDw,dRdDw)
end
type InertialPose3Container <: PreintContainer
  rRp::Array{Float64,2}
  rPosp::Vector{Float64}
  rVelp::Vector{Float64}
  pBw::Vector{Float64}
  pBa::Vector{Float64}
  rnTime::Int64
  InertialPose3Container(;
        rRp::Array{Float64,2}=eye(3),
        rPosp::Vector{Float64}=zeros(3),
        rVelp::Vector{Float64}=zeros(3),
        pBw::Vector{Float64}=zeros(3),
        pBa::Vector{Float64}=zeros(3),
        rnTime::Int64=0  ) = new( rRp,rPosp,rVelp,pBw,pBa,rnTime )
end

function vee(ip3::InertialPose3Container)
  ret = zeros(15)
  ret[1:6] = veeEuler(  SE3(ip3.rPosp, SO3(rRp))  )

end

function oplus(xi::InertialPose3Container, Dx::InertialPose3Container)
  return InertialPose3Container(
        xi.rRp*Dx.rRp,
        xi.rPosp + xi.rRp*Dx.rPosp,
        xi.rVelp + xi.rRp*Dx.rVelp*(Dx.rnTime*1e-9),
        xi.pBw + Dx.pBw,
        xi.pBa + Dx.pBa,
        xi.rnTime+Dx.rnTime  )
end
⊕(xi::InertialPose3Container, Dx::InertialPose3Container) = oplus(xi,Dx)

function zetaEmbedding(posei::InertialPose3Container, posej::InertialPose3Container; rGrav=[0.0;0.0;9.81])
  z = zeros(30)
  # z[1:3] = logmap(SO3(posei.rRp'*posej.rRp))

  z[4:6] = posej.pBw
  z[7:9] = posej.rVelp
  z[10:12] = posej.rPosp
  z[13:15] = posej.pBa
  z[16:18] = posei.pBw
  z[19:21] = posei.rVelp
  z[22:24] = posei.rPosp
  z[25:27] = posei.pBa
  z[28:30] = rGrav
  return z
end

# predict the preintegral value
function constructL(posei::InertialPose3Container, Dt::Float64)
  biRw = posei.rRp'
  L = zeros(15,30)
  L[1:3,1:3] = eye(3)
  L[7:9,7:9] = biRw
  L[10:12,10:12] = biRw
  L[7:9,19:21] = -biRw
  L[10:12,19:21] = -biRw*Dt # preintegral constant velocity component in delta position
  L[10:12,22:24] = -biRw
  return L
end

# First term in Taylor expansion
function constructC1{T <: PreintContainer}(posei::InertialPose3Container, pido::T, Dt::Float64)
  biRw = posei.rRp'
  C1 = zeros(15,30)
  g1 = -biRw*Dt
  g2 = 0.5*g1*Dt
  C1[4:6,4:6] = eye(3)
  C1[13:15,13:15] = eye(3)
  C1[4:6,16:18] = -eye(3)
  C1[13:15,25:27] = -eye(3)
  C1[7:9,28:30] = g1
  C1[10:12,28:30] = g2
  C1[1:3,16:18] = pido.dRdDw
  C1[7:9,16:18] = pido.dVdDw
  C1[10:12,16:18] = pido.dPdDw
  C1[7:9,25:27] = pido.dVdDa
  C1[10:12,25:27] = pido.dPdDa
  return C1
end

function preintMeas{T <: PreintContainer}(pido::T)
  return [logmap(SO3(pido.iRj));zeros(3);pido.iDvj;pido.iDppj;zeros(3)] # temporarily suppressing bias updates
end

function predictDeltaXij{T <: PreintContainer}(pido::T, posei::InertialPose3Container, posej::InertialPose3Container; rGrav=[0.0;0.0;9.81])
  zet = zetaEmbedding(posei, posej, rGrav=rGrav)
  Dt = Float64(posej.rnTime - posei.rnTime)*1e-9
  L = constructL(posei,Dt)
  C1 = constructC1(posei, pido, Dt)
  return (L-C1)*zet #-0.5*C2*zet.^2
end

function residual!(res::Vector{Float64},
        pioc::InertialPose3Container,
        picg::PreintegralCompensationGradients,
        posei::InertialPose3Container,
        posej::InertialPose3Container, rGrav=[0.0;0.0;9.81]  )
  #
  res[1:15] = preintMeas(pioc) - predictDeltaXij(picg, posei, posej, rGrav=rGrav)
  nothing
end


"""
    InertialPose3(zij::Distributions.MvNormal, pioc::InertialPose3Container, pido::PreintegralCompensationGradients)

Zij is entropy of veeLie15, pioc is preintegral measurements, pido is compensation gradients.
"""
type InertialPose3 <: RoME.BetweenPoses
  Zij::Distribution
  pioc::InertialPose3Container
  picg::PreintegralCompensationGradients
  reuse::Vector{Tuple{InertialPose3Container,InertialPose3Container, InertialPose3Container}} # number of threads

  InertialPose3(
    zij::Distributions.MvNormal,
    pioc::InertialPose3Container,
    picg::PreintegralCompensationGradients  ) =
      new( zij, pioc, picg,
        fill( (InertialPose3Container(),
               InertialPose3Container(),
               InertialPose3Container()),
        Threads.nthreads() )
      )
end

function getSample(ip3::InertialPose3, N::Int=1)
  return (rand( ip3.Zij, N )', )
end

function (ip3::InertialPose3)(
        res::Vector{Float64},
        idx::Int,
        meas::Tuple,
        wIPi::Array{Float64,2},
        wIPj::Array{Float64,2}  )
  #
  # Function can be massively improved. Just getting it all wired at first.
  # get pointer to memory, local to this thread
  posei, posej, ENT = ip3.reuse[Threads.threadid()]

  ## this part must be moved in the general ApproxConv.jl code, since this reassignment
  ## does not have to happen every time
  # repoint to existing values for first pose
  posei.rRp[1:3,1:3] = convert(SO3, Euler(wIPi[4:6, idx]...)).R
  posei.rPosp = wIPi[1:3, idx]
  posei.rVelp = wIPi[7:9, idx]
  posei.pBw = wIPi[10:12, idx]
  posei.pBa = wIPi[13:15, idx]

  # repoint to existing load values for second pose
  posei.rRp[1:3,1:3] = convert(SO3, Euler(wIPj[4:6, idx]...)).R
  posei.rPosp = wIPj[1:3, idx]
  posei.rVelp = wIPj[7:9, idx]
  posei.pBw = wIPj[10:12, idx]
  posei.pBa = wIPj[13:15, idx]

  # repoint to existing load values for second pose
  ENT.rRp[1:3,1:3] = convert(SO3, so3(meas[1][4:6, idx])).R,
  ENT.rPosp = meas[1][1:3, idx],
  ENT.rVelp = meas[1][7:9, idx],
  ENT.pBw = meas[1][10:12, idx],
  ENT.pBa = meas[1][13:15, idx],
  ENT.rnTime = ip3.pioc.rnTime # need time to enact velocity noise

  keeprntime = deepcopy(posej.rnTime)
  noisyposej = posej⊕ENT
  poisyposej.rnTime =keeprntime  # rest timestamp to correct value
  residual!(res, ip3.picg, posei,  noisyposej)
  nothing
end

type PackedInertialPose3 <: IncrementalInference.PackedInferenceType
  vecZij::Array{Float64,1} # 3translations, 3rotation, 3 velocities
  vecCov::Array{Float64,1}
  dimc::Int64
  vecpioc::Vector{Float64}
  rnTime::Int64
  picgvecdPdDa::Vector{Float64}
  picgvecdVdDa::Vector{Float64}
  picgvecdPdDw::Vector{Float64}
  picgvecdVdDw::Vector{Float64}
  picgvecdRdDw::Vector{Float64}
  PackedInertialPose3() = new()
  PackedInertialPose3(ip3::InertialPose3) = new( ip3.Zij.μ, ip3.Zij.Σ.mat[:], size(ip3.Zij.Σ.mat,1),
            veeQuaternion(ip3.pioc), ip3.pioc.rnTime,
            ip3.dPdDa[:] ,ip3.dVdDa[:] ,ip3.dPdDw[:] ,ip3.dVdDw[:] ,ip3.dRdDw[:] )
end

convert(::Type{PackedInertialPose3}, ip3::InertialPose3) = PackedInertialPose3(ip3)

function convert(::Type{InertialPose3}, pip3::PackedInertialPose3)
  pioc = InertialPose3Container(
      rPosp=pip3.vecpioc[1:3],
      rRp=Quaternion(vecpioc[4],vecpioc[5:7]),
      rVelp=vecpioc[8:10],
      pBw=vecpioc[11:13],
      pBa=vecpioc[14:16],
      rnTime=pip3.rnTime
  )
  picg = PreintegralCompensationGradients(
      reshapeVec2Mat(pip3.picgvecdPdDa, 3),
      reshapeVec2Mat(pip3.picgvecdVdDa, 3),
      reshapeVec2Mat(pip3.picgvecdPdDw, 3),
      reshapeVec2Mat(pip3.picgvecdVdDw, 3),
      reshapeVec2Mat(pip3.picgvecdRdDw, 3)
  )
  InertialPose3(Distributions.MvNormal(pip3.vecZij, reshapeVec2Mat(pip3.vecCov, pip3.dimc)), pioc, picg)
end




function compare(a::PreintegralCompensationGradients, b::PreintegralCompensationGradients)
  TP = true
  TP = TP && norm( a.dPdDa - b.dPdDa ) < 1e-10
  TP = TP && norm( a.dVdDa - b.dVdDa ) < 1e-10
  TP = TP && norm( a.dPdDw - b.dPdDw ) < 1e-10
  TP = TP && norm( a.dVdDw - b.dVdDw ) < 1e-10
  TP = TP && norm( a.dRdDw - b.dRdDw ) < 1e-10
  return TP
end
function compare(a::InertialPose3, b::InertialPose3)
  TP = true
  TP = TP && norm(a.Zij.μ - b.Zij.μ) < 1e-10
  TP = TP && norm(a.Zij.Σ.mat - b.Zij.Σ.mat) < 1e-10
  TP = TP && compare(a.pido, b.pido)
  return TP
end









# we also need a prior for InertialPose3



type PriorInertialPose3 <: IncrementalInference.FunctorSingleton
  Zi::Distribution
end
function getSample(prip3::PriorInertialPose3, N::Int=1)
  return (rand( prip3.Zi, N )', )
end


type PackedPriorInertialPose3 <: IncrementalInference.PackedInferenceType
  vecZi::Array{Float64,1} # 3translations, 3rotation, 3 velocities
  vecCov::Array{Float64,1}
  dimc::Int64
  PackedPriorInertialPose3() = new()
  PackedPriorInertialPose3(prip3::PriorInertialPose3) = new( prip3.Zi.μ, prip3.Zi.Σ.mat[:], size(prip3.Zi.Σ.mat,1)  )
end

convert(::Type{PackedPriorInertialPose3}, prip3::PriorInertialPose3) = PackedPriorInertialPose3(prip3)

function convert(::Type{PriorInertialPose3}, pprip3::PackedPriorInertialPose3)
  PriorInertialPose3(Distributions.MvNormal(pprip3.vecZi, reshapeVec2Mat(pprip3.vecCov, pprip3.dimc)) )
end




#
