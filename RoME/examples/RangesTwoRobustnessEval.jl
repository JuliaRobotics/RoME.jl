# draw robustness evaluation pictures

using DataFrames, Gadfly


folderloc = "/home/dehann/Videos/"

NN = Int[25,50,75,100,125,150,175,200,250,300,350,400,500,600,700,800,900,1000]

NN = Int[25,50,150,350,700]

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
particles = N*ones(length(df[:x1].data[3:end])),
Particles="$(N)"
))



nPE = "slamedonutAllPE_$(N).txt"
fi = joinpath(folderloc, nPE)
df = readtable(fi, header=false)

push!(allerrors, DataFrame(
variable = union(["L$(i)" for i in 3:4],["X$(i)" for i in 1:13]),
distance = df[:x1].data[3:end],
particles = N*ones(length(df[:x1].data[3:end])),
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

vDF = vcat(alllikelihoods[end:-1:1]...)

pl = plot(vDF,Geom.point,color=:Particles,
x=:variable,y=:likelihood,
Guide.xticks(orientation=:vertical))

draw(PDF("slamedonutLikelihoods.pdf",12cm,7cm),pl)

vDF = vcat(alllikelihoods...)

pl = plot(vDF, Geom.boxplot,
x=:particles, y=:likelihood,
Coord.cartesian(xmax=1025))

draw(PDF("slamedonutallLK.pdf",12cm,7cm),pl)


vDF = vcat(allerrors[end:-1:1]...)

pl = plot(vDF,Geom.point,color=:Particles,
x=:variable,y=:distance,
Guide.xticks(orientation=:vertical))

draw(PDF("slamedonutPoseErrs.pdf",12cm,7cm),pl)

vDF = vcat(allerrors...)

pl = plot(vDF, Geom.boxplot,
x=:particles, y=:distance,
Coord.cartesian(xmax=1025))

draw(PDF("slamedonutallPE.pdf",12cm,7cm),pl)


DFcomp = DataFrame(
particles=NN,
seconds=computetime,
time="measured"
)
DFcompE = DataFrame(
particles=NN,
seconds=computetime/5,
time="possible"
)

DFc = vcat(DFcomp, DFcompE)


pl = plot(DFc,x=:particles,y=:seconds, color=:time, Geom.line)

draw(PDF("slamedonutcomputetimes.pdf",12cm,7cm),pl)







#
