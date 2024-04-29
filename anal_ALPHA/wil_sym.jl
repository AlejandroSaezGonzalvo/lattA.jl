#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "N300"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

pp_sym, ap_sym, corr, corr_val, corrw, dSdm, w = read_ens_wil(path, ens, legacy=true)

#======== compute observables ========#

if ens.id == "N300"
    mpi = get_m_ALPHA(pp_sym[1], ens, "pion_wil", pl=true, tm=35, wpm=wpm)
elseif ens.id == "N202"
    mpi = get_m_ALPHA(pp_sym[1], ens, "pion_wil", pl=true, tm=32, wpm=wpm)
else
    mpi = get_m_ALPHA(pp_sym[1], ens, "pion_wil", pl=true, wpm=wpm)
end
mk = mpi
m12 = get_mpcac_ALPHA(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, wpm=wpm)
m13 = m12
fpi = get_f_wil_ALPHA(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, wpm=wpm)
fk = fpi

ZA = beta_ZA[ens.beta]
bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = ZA * (1 + bAtil * m12) * fpi
fk = ZA * (1 + bAtil * m13) * fk

m12_I = (1 + beta_bap[ens.beta] * m12) * m12

#======== compute t0/a² ===============#

if ens.id == "H400"
    t0, YW, WY = get_t0_ALPHA(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=10, pl=true)
elseif ens.id == "N202"
    t0, YW, WY = get_t0_ALPHA(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=29, pl=true)
elseif ens.id == "N300"
    t0, YW, WY = get_t0_ALPHA(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=33, pl=true)
else
    t0, YW, WY = get_t0_ALPHA(path, ens, [40,60], rw=true, info=true, wpm=wpm, pl=true)
end

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_obs_wil_un_ALPHA.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

#============ get md & dm =============#

phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
a = phi4
md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
v1 = v2 = s1 = s2 = 0
for i in 1:length(md_v)
    v1 += md_v[i][1]
    v2 += md_v[i][2]
    s1 += md_s[i][1]
    s2 += md_s[i][2]
end
s1 += md_s[end][1]
s2 += md_s[end][2]
dm = (phi4_ph - phi4) / (2*s1 + s2 + v1 + v2)

fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/dm_", ens.id, "_phi4=", round(value(phi4_ph), digits=5), "_ALPHA.bdio"), "w")
write_uwreal(dm, fb, 1)
BDIO_close!(fb)

if md_meas == true
    obs_md = Array{uwreal,1}()
    for a in [phi4; obs]
        md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
        md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
        v1 = v2 = s1 = s2 = 0
        for i in 1:length(md_v)
            v1 += md_v[i][1]
            v2 += md_v[i][2]
            s1 += md_s[i][1]
            s2 += md_s[i][2]
        end
        s1 += md_s[end][1]
        s2 += md_s[end][2]
        push!(obs_md, 2*s1 + s2 + v1 + v2)
    end

    ## now compute only strange derivatives wrt phi4 to interpolate as in Ben Strassberger's Thesis
    phi2 = 8 * t0 * mpi ^ 2
    t0fpik = sqrt(8 * t0) * 2/3 * (fk + 0.5 * fpi)
    for a in [phi4; phi2; t0fpik]
        md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
        md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
        v2 = s2 = 0
        for i in 1:length(md_v)
            v2 += md_v[i][2]
            s2 += md_s[i][2]
        end
        s2 += md_s[end][2]
        push!(obs_md, s2 + v2)
    end
    obs_md[end-1] = obs_md[end-1] / obs_md[end-2]; obs_md[end] = obs_md[end] / obs_md[end-2]

    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_md_wil_ALPHA.bdio"), "w")
    for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
    BDIO_close!(fb)
end
