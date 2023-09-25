import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#======== read correlators ===========#

pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=true);
pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(pp)-1];
ap, apw, w = get_corr_wil(path, ens, "G5", "G0G5", rw=true, info=true, legacy=true);
ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)-1];

#======== compute observables ========#

mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true)
mk = mpi
m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true)
m13 = m12
fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true)
fk = fpi

mpi, mk, m12, m13, fpi, fk = mpi[1], mk[1], m12[1], m13[1], fpi[1], fk[1]
mpi, fpi, fk = fve(mpi, mk, fpi, fk, ens)

#======== compute t0/a² ===============#

t0 = get_t0(path, ens, [40,60], rw=true, wpm=wpm)

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)







