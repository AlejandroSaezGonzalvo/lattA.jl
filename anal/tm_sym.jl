import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

#include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

pp_sym, ap_sym, corrw, dSdm = read_ens_tm_sym(path, ens, legacy=true)

#======== compute observables ========#

tm = [collect(10:5:div(ens.T,2)-5), collect(10:5:div(ens.T,2)-5)]
tM = [collect(div(ens.T,2)+5:5:ens.T-10), collect(div(ens.T,2)+5:5:ens.T-10)]

mpi = Array{uwreal,1}()
m12 = Array{uwreal,1}()
fpi = Array{uwreal,1}()
for i in 1:length(pp_sym)
    mpi_aux = get_m(pp_sym[i], ens, "pion_tm", tm=tm, tM=tM)
    push!(mpi, mpi_aux[1])
    m12_aux = get_mpcac(pp_sym[i], ap_sym[i], ens, "pion_tm", tm=tm, tM=tM)
    push!(m12, m12_aux[1])
    fpi_aux = get_f_tm(pp_sym[i], mpi[i], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
    push!(fpi, fpi_aux[1])
end
fk = fpi
mk = mpi

for i in 1:length(mpi) mpi[i], fpi[i], fk[i] = fve(mpi[i], mk[i], fpi[i], fk[i], ens) end
fk = fpi ## need to "impose" this after fve in case of sym ens
mk = mpi

#======== mass shift with fit =========#

par = Array{uwreal,1}()
fb = BDIO_open("/home/asaez/cls_ens/results/der_1q.bdio", "r")
BDIO_seek!(fb); push!(par, read_uwreal(fb))
while BDIO_seek!(fb, 2) == true push!(par, read_uwreal(fb)) end 
BDIO_close!(fb)

obs = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "r")
BDIO_seek!(fb); push!(obs, read_uwreal(fb))
for i in 2:3 BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb)) end
t0, mpi_w, mk_w = obs
phi4_w = 8 * t0 * (mk_w ^ 2 + 0.5 * mpi_w ^ 2)
phi2_w = 8 * t0 * mpi_w ^ 2

ZP = beta_ZP[ens.beta]
phi2 = [8 * t0 * mpi[i] ^ 2 for i in 1:length(mpi)]
phi4 = [8 * t0 * (mk[i] ^ 2 + 0.5 * mpi[i] ^ 2) for i in 1:length(mpi)]
fpik = [sqrt(t0) * 2/3 * (fk[i] + 0.5 * fpi[i]) for i in 1:length(mpi)]

fpik_sh = [fpik[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[10:12]) for i in 1:length(fpik)]
phi2_sh = [phi2[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[13:15]) for i in 1:length(phi2)]
phi4_sh = [phi4[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[16:18]) for i in 1:length(phi4)]
m12_sh = [m12[i] + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[19:23]) * ZP / sqrt(t0) for i in 1:length(m12)]

#======== save BDIO ===================#

obs = [mpi, m12, fpi]
obs_str = ["mpi", "m12", "fpi"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "w")
    for i in 1:length(obs[j]) write_uwreal(obs[j][i], fb, i) end
    BDIO_close!(fb)
end

obs = [phi2_sh, phi4_sh, m12_sh, fpik_sh]
obs_str = ["phi2", "phi4", "m12", "fpik"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_", obs_str[j], "_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
    for i in 1:length(obs[j]) write_uwreal(obs[j][i], fb, i) end
    BDIO_close!(fb)
end

#============ get md ==================#

if md_meas == true
    obs_md = [Array{uwreal,1}() for a in obs]
    for i in 1:length(obs)
        for a in obs[i]
            md_s = [md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]
            s1 = s2 = 0
            for i in 1:length(md_s)
                s1 += md_s[i][1]
                s2 += md_s[i][2]
            end
            push!(obs_md[i], 2*s1 + s2)
        end
    end


    for j in 1:length(obs_st)
        fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_", obs_st[j],"_md_tm.bdio"), "w")
        for i in 1:length(obs[j]) write_uwreal(obs_md[j][i], fb, i) end
        BDIO_close!(fb)
    end
end







