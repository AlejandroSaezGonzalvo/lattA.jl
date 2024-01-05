import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H102r001"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=true);
pp_sym = [corr_sym(pp[i], pp[i+1], +1) for i in 1:2:length(pp)-1];
ap, apw, w = get_corr_wil(path, ens, "G5", "G0G5", rw=true, info=true, legacy=true);
ap_sym = [corr_sym(ap[i], ap[i+1], -1) for i in 1:2:length(ap)-1];

pp_d1 = get_corr_wil(path, ens, "G5_d1", "G5_d1", rw=true, legacy=true);
pp_d2 = get_corr_wil(path, ens, "G5_d2", "G5_d2", rw=true, legacy=true);
ap_d1 = get_corr_wil(path, ens, "G5_d1", "G0G5_d1", rw=true, legacy=true);
ap_d2 = get_corr_wil(path, ens, "G5_d2", "G0G5_d2", rw=true, legacy=true);
dSdm = get_dSdm(path, ens)

pp_val = [[pp_d1[i], pp_d2[i]] for i in 1:length(pp_d1)];
ap_val = [[ap_d1[i], ap_d2[i]] for i in 1:length(ap_d1)];
corr = [[pp[i] for i in 1:length(pp)]; [ap[i] for i in 1:length(ap)]];
corr_val = [[pp_val[i] for i in 1:length(pp)]; [ap_val[i] for i in 1:length(ap)]];
corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

#======== compute observables ========#

tm = [collect(10:5:div(ens.T,2)-5), collect(10:5:div(ens.T,2)-5)]
tM = [collect(div(ens.T,2)+5:5:ens.T-10), collect(div(ens.T,2)+5:5:ens.T-10)]

mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
mk = get_m(pp_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=true, wpm=wpm, tm=tm, tM=tM)

mpi, mk, m12, m13, fpi, fk = mpi[1], mk[1], m12[1], m13[1], fpi[1], fk[1]
mpi, fpi, fk = fve(mpi, mk, fpi, fk, ens)

bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = (1 + bAtil * m12) * fpi
fk = (1 + bAtil * m13) * fk

#======== compute t0/a² ===============#

t0, YW, WY = get_t0(path, ens, [40,80], rw=true, info=true, wpm=wpm)

#======== mass shift with fit =========#

par = Array{uwreal,1}()
fb = BDIO_open("/home/asaez/cls_ens/results/der_1q.bdio", "r")
BDIO_seek!(fb); push!(par, read_uwreal(fb))
while BDIO_seek!(fb, 2) == true push!(par, read_uwreal(fb)) end 
BDIO_close!(fb)

ZA = beta_ZA[ens.beta]
phi2 = 8 * t0 * mpi ^ 2
phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
fpik = ZA * sqrt(t0) * 2/3 * (fk + 0.5 * fpi)

fpik_sh = fpik + (phi4_ph - phi4) * md_1(phi2, 1/t0, par[1:3])
phi2_sh = phi2 + (phi4_ph - phi4) * md_1(phi2, 1/t0, par[4:6])
t0_sh = t0 + (phi4_ph - phi4) * md_1(phi2, 1/t0, par[7:9])

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

obs = [t0_sh, phi2_sh, fpik_sh]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

#============ get md ==================#

if md_meas == true
    phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
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
        push!(obs_md, s2 + v2)
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







