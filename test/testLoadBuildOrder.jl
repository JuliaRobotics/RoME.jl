# test load order

using Caesar


fg31bd = initfg()
loadDFG(joinpath(ENV["HOME"],"Downloads/test/fg_beforedownsolve"), Main, fg31bd)

sfg31bd = buildSubgraphFromLabels(fg31bd, [:x167;:x168;:x169;:x170])

fct = getFactor(sfg31bd, :x168x169f1)
@show getFactorType(fct)
getData(fct).fncargvID
