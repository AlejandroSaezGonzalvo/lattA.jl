import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H400"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

pp_sym, ap_sym, corr, corr_val, corrw, dSdm, w = read_ens_wil(path, ens, legacy=true)

#======== compute observables ========#

tm = [[10], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
tM = [[ens.T-10], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]

mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, tm=tm, tM=tM, wpm=wpm)
mk = mpi
m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, tm=tm, tM=tM, wpm=wpm)
m13 = m12
fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, tm=tm, tM=tM, wpm=wpm)
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

fb = BDIO_open(string("/home/asaez/cls_ens/results/shifted/dm_", ens.id, "_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
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

    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_md_wil.bdio"), "w")
    for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
    BDIO_close!(fb)
end
