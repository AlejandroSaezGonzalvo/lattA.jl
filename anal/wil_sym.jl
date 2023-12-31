import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

pp_sym, ap_sym, corr, corr_val, corrw, dSdm = read_ens_wil(path, ens, legacy=true)

#======== compute observables ========#

tm = [collect(10:5:div(ens.T,2)-5), collect(10:5:div(ens.T,2)-5)]
tM = [collect(div(ens.T,2)+5:5:ens.T-10), collect(div(ens.T,2)+5:5:ens.T-10)]

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

bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = (1 + bAtil * m12) * fpi
fk = (1 + bAtil * m13) * fk

#======== compute t0/a² ===============#

t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm)

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
phi2_sh = phi2 + (phi4_ph - phi4) * 2 / 3 ## symmetric point definition
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
