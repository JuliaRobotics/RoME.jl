# draw robustness evaluation pictures

using DataFrames, Gadfly


folderloc = "/home/dehann/Videos/"

NN = Int64[25,50,75,100,125,150,175,200,250,300,350,400]

NN = Int64[25,50,200,400]

# N = NN[12]


alllikelihoods = DataFrame[]
allerrors = DataFrame[]

computetime = Float64[]

for N in NN


nLK = "slamedonutAllLK_$(N).txt"
fi = joinpath(folderloc, nLK)
df = readtable(fi, header=false)


push!(alllikelihoods, DataFrame(
variable = union(["L$(i)" for i in 3:4],["X$(i)" for i in 1:13]),
likelihood = df[:x1].data[3:end],
Particles="$(N)"
))



nPE = "slamedonutAllPE_$(N).txt"
fi = joinpath(folderloc, nPE)
df = readtable(fi, header=false)

push!(allerrors, DataFrame(
variable = union(["L$(i)" for i in 3:4],["X$(i)" for i in 1:13]),
distance = df[:x1].data[3:end],
Particles="$(N)"
))



nPE = "slamedonutBatchTime_$(N).txt"
fi = joinpath(folderloc, nPE)
fid = open(fi,"r")
vals = readline(fid)
close(fid)
val = parse(Float64,vals)
push!(computetime, val)

end

vDF = vcat(alllikelihoods...)

plot(vDF,Geom.point,color=:Particles,
x=:variable,y=:likelihood)



vDF = vcat(allerrors...)

plot(vDF,Geom.point,color=:Particles,
x=:variable,y=:distance)



DFcomp = DataFrame(
particles=NN,
seconds=computetime,
Time="measured"
)
DFcompE = DataFrame(
particles=NN,
seconds=computetime/5,
Time="expected"
)

DFc = vcat(DFcomp, DFcompE)



plot(DFc,x=:particles,y=:seconds, color=:Time, Geom.line)



#
