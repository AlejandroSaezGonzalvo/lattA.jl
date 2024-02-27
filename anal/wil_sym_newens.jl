import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "J500"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

if id == "N302"
    pp_sym, ap_sym = read_ens_TSM(path, ens) 
else
    pp_sym, ap_sym = get_corr_TSM_multichunks(path, ens)
end

#======== compute observables ========#

tm = [[10], collect(20:10:div(ens.T,2)-10)]
tM = [[ens.T-10], collect(div(ens.T,2)+10:10:ens.T-20)]

mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, tm=tm, tM=tM)
mk = mpi
m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, tm=tm, tM=tM)
m13 = m12
fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, tm=tm, tM=tM)
fk = fpi

mpi, mk, m12, m13, fpi, fk = mpi[1], mk[1], m12[1], m13[1], fpi[1], fk[1]
mpi, fpi, fk = fve(mpi, mk, fpi, fk, ens)
fk = fpi ## need to "impose" this after fve in case of sym ens
mk = mpi

ZA = beta_ZA[ens.beta]
bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = ZA * (1 + bAtil * m12) * fpi
fk = ZA * (1 + bAtil * m13) * fk

#======== compute t0/a² ===============#

t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM, pl=true)

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)