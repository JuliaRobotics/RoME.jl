# keep data from Preintegration process

using TransformUtils
using PyCall
using PyLCM

using RoME

function gen_bindings()
    @show lcmtpath = joinpath(dirname(@__FILE__),"lcmtypes")
    run(`lcm-gen -p --ppath $(lcmtpath) $(lcmtpath)/devjl_preint_d33_short_t.lcm`)
    println("Adding lcmtypes dir to Python path: $(lcmtpath)")
    unshift!(PyVector(pyimport("sys")["path"]),lcmtpath)
end

gen_bindings();

@pyimport devjl

function handlepreints(channel, msgdata, storedata)
  msg = devjl.preint_d33_short_t[:decode](msgdata)
  @show msg[:utime]
  @show iDtj = msg[:Dt]
  @show iDppj= Float64[ msg[:iDppj]... ]
  @show iDvj = Float64[ msg[:iDvj]... ]
  iRj = reshape( Float64[ msg[:iDRj]... ], 3,3 )
  @show iDqj = convert(Quaternion, SO3(iRj))

  msgCov = reshape( Float64[ msg[:Cov]... ], 16,16 )

  P = zeros(15,15)
  P[1:12,1:12]    = msgCov[1:12,1:12]
  P[1:12,13:15]   = msgCov[1:12,14:16]
  P[13:15,1:12]   = msgCov[14:16,1:12]
  P[13:15,13:15]  = msgCov[14:16,14:16]
  @show dDw = Float64[ msg[:dDw]... ]
  @show dDa = Float64[ msg[:dDa]... ]

  pioc = InertialPose3Container(
      rRp=iRj,
      rPosp=iDppj,
      rVelp=iDvj,
      pBw=dDw,
      pBa=dDa,
      rnTime=round(Int, iDtj*1000000000)
  )

  picg = PreintegralCompensationGradients(
      dPdDa = reshape( Float64[msg[:dPdDa]...] , 3,3 ),
      dVdDa = reshape( Float64[msg[:dVdDa]...] , 3,3 ),
      dPdDw = reshape( Float64[msg[:dPdDw]...] , 3,3 ),
      dVdDw = reshape( Float64[msg[:dVdDw]...] , 3,3 ),
      dRdDw = reshape( Float64[msg[:dRdDw]...] , 3,3 )
  )

  push!(storedata, (pioc, picg, P) )

  nothing
end


lc = LCM()

DATA = Vector{Tuple{InertialPose3Container,PreintegralCompensationGradients, Array{Float64,2}}}()

store = (c, m) -> handlepreints(c, m, DATA)

subscribe(lc, "IMU_PREINTEGRALS", store)

iters=9
for i in 1:iters
# while true
  handle(lc)
end


using JLD

file = joinpath(dirname(@__FILE__),"preintstationarydata.jld")
@save file DATA
# Base.rm(file)





#
