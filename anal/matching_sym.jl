#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#=========== read bdio ================#

obs = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ens.id, "_obs_wil_un.bdio"), "r")
BDIO_seek!(fb); push!(obs, read_uwreal(fb))
for i in 2:7 BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb)) end
BDIO_close!(fb)
t0, mpi_w, mk_w, m12_w, m13_w, fpi_w, fk_w = obs
phi4_w = 8 * t0 * (mk_w ^ 2 + 0.5 * mpi_w ^ 2)
phi2_w = 8 * t0 * mpi_w ^ 2
uwerr(phi2_w)
uwerr(phi4_w)

obs = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
obs_str = ["mpi", "m12", "fpi"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "r")
    BDIO_seek!(fb); push!(obs[j], read_uwreal(fb))
    while BDIO_seek!(fb, 2) == true push!(obs[j], read_uwreal(fb)) end 
    BDIO_close!(fb)
end
mpi = obs[1]
m12 = obs[2] ## not improved, not renormalized
fpi = obs[3] ## not improved, not renormalized

ZA = beta_ZA[ens.beta]
bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = ZA * (1 + bAtil * m12) * fpi

phi2 = [8 * t0 * mpi[i] ^ 2 for i in 1:length(mpi)]
phi4 = [8 * t0 * (mpi[i] ^ 2 + 0.5 * mpi[i] ^ 2) for i in 1:length(mpi)]
fpik = [sqrt(t0) * 2/3 * (fpi[i] + 0.5 * fpi[i]) for i in 1:length(mpi)]
uwerr.(phi2)
uwerr.(phi4)
uwerr.(m12)
uwerr.(fpik)

mpi_w, mk_w, fpi_w, fk_w = fve(mpi_w, mk_w, fpi_w, fk_w, ens)
for i in 1:9 mpi[i], a, fpi[i], b = fve(mpi[i], mpi[i], fpi[i], fpi[i], ens) end

#========== mass shift =============#

par = Array{uwreal,1}()
fb = BDIO_open("/home/asaez/cls_ens/results/derivatives/der_1q.bdio", "r")
BDIO_seek!(fb); push!(par, read_uwreal(fb))
while BDIO_seek!(fb, 2) == true push!(par, read_uwreal(fb)) end 
BDIO_close!(fb)

ZA = beta_ZA[ens.beta]
ZP = beta_ZP[ens.beta]

phi4_w = 8 * t0 * (mk_w ^ 2 + 0.5 * mpi_w ^ 2)
phi2_w = 8 * t0 * mpi_w ^ 2
uwerr(phi2_w)
uwerr(phi4_w)
fpik_w = sqrt(t0) * 2/3 * (fk_w + 0.5 * fpi_w)

fpik_w_sh = fpik_w + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[1:3]) ## already improved and renormalized
phi2_w_sh = phi2_w + (phi4_ph - phi4_w) * 2 / 3 ## symmetric point definition
t0_sh = t0 + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[7:9])
m12_w_sh = m12_w + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[30:34]) * ZP / sqrt(t0) / ZA ## not improved, not renormalized
m12_w_sh_I = (1 + beta_bap[ens.beta] * m12_w_sh) * m12_w_sh
m12_w_I = (1 + beta_bap[ens.beta] * m12_w) * m12_w
fpi_w_sh = fpi_w + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[24:26]) / sqrt(t0)

fpik_sh = [fpik[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[10:12]) for i in 1:length(fpik)]
fpi_sh = [fpi[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[40:42]) for i in 1:length(fpi)]
phi2_sh = [phi2[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[13:15]) for i in 1:length(phi2)]
phi4_sh = [phi4[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[16:18]) for i in 1:length(phi4)]
m12_sh = [m12[i] + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[19:23]) * ZP / sqrt(t0) / ZA for i in 1:length(m12)]
uwerr.(fpik_sh)
uwerr.(fpi_sh)
uwerr.(phi2_sh)
uwerr.(phi4_sh)
uwerr.(m12_sh)

#======== match & full twist =======#

y = [m12_sh, 2/3 .* phi4_sh]
target = [0 .* m12_sh, [2/3 * phi4_ph for i in 1:length(phi4)]]
y = y .- target

kappa, mul = ens_kappa[id], ens_mul[id]
x = [[[kappa[1] for i in 1:3]; [kappa[2] for i in 1:3]; [kappa[3] for i in 1:3]] [mul; mul; mul]]

up, chi2, chi_exp, pv = fit_alg_LBFGS([match_m12_sym, match_phi2_sym],[x, x],y,6,[kappa[2], mul[2]],wpm=wpm) ##kappa->up[1], mul->up[2]

matching_sym_plot()

#========= interpolate fpik ========#

y = fpik_sh
cc = 0
while cc == 0
    try
        global up_fpik, chi2, chi_exp, pv = fit_alg(interp_fpik_sym,x,y,4,ens_up_fpik[ens.id] .+ rand(5),wpm=wpm)
        global fpik_matched = interp_fpik_sym([up[1] up[2]],up_fpik)[1]
        uwerr(fpik_matched)
        cc+=1
    catch e
    end
end
interp_fpik_sym_plot()

#========= save bdio ===============#

obs = [t0_sh, phi2_w_sh, m12_w_sh_I, m12_w_sh_I, fpi_w_sh, fpi_w_sh, fpik_w_sh]
fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted/", ens.id, "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

obs = [up[1], up[2], up[2], fpik_matched / sqrt(t0_sh), fpik_matched / sqrt(t0_sh), fpik_matched]
fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted/", ens.id, "_obs_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)