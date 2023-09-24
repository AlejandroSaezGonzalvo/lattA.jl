import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
wpm = Dict{String, Vector{Float64}}()
wpm[id] = [-1.0, -1.0, 4.0, -1.0]
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=true);
pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(pp)-1];
ap, apw, w = get_corr_wil(path, ens, "G5", "G0G5", rw=true, info=true, legacy=true);
ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)-1];

mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true)
mk = mpi
println("mpi = ", mpi[1])
m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true)
m13 = m12
println("m12 = ", m12[1])
fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true)
fk = fpi
println("fpi = ", fpi[1])

#mpi, fpi, fk = fve(mpi[1], mk[1], fpi[1], fk[1], ens)

t0, YW, WY = get_t0(path, ens, [40,60], dtr=2, rw=true, info=true)







