import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#======== read correlators ===========#

pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=true);
pp_sym = [corr_sym(pp[i], pp[i+9], +1) for i in 1:2:length(pp)-1];
ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=true);
ap_sym = [corr_sym(ap[i], ap[i+9], -1) for i in 1:2:length(ap)-1];

dSdm = get_dSdm(path, ens)

corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

#======== compute observables ========#

mpi = fpi = m12 = Array{uwreal,1}()
for i in 1:length(pp_sym)
    mpi_aux = get_m(pp_sym[1], ens, "pion_tm")
    push!(mpi, mpi_aux[1])
    m12_aux = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_tm")
    push!(m12, m12_aux[1])
    fpi_aux = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_tm")
    push!(fpi, fpi_aux[1])
end

mpi, fpi, fk = fve.(mpi, mk, fpi, fk, ens)
fk = fpi ## need to "impose" this after fve in case of sym ens

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_tm_un.bdio"), "w")
for a in obs write_uwreal(a, fb, i) end
BDIO_close!(fb)

#============ get md ==================#

obs_md = Array{uwreal,1}()
for a in obs
    md_s = [md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]
    s1 = s2 = 0
    for i in 1:length(md_v)
        s1 += md_s[i][1]
        s2 += md_s[i][2]
    end
    push!(obs_md, 2*s1 + s2)
end

#======== save BDIO ===================#

fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_md_tm.bdio"), "w")
for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
BDIO_close!(fb)







