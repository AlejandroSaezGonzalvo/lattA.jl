#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "H105r005"
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

obs = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
obs_str = ["mpi", "mk", "m12", "m13", "fpi", "fk"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "r")
    BDIO_seek!(fb); push!(obs[j], read_uwreal(fb))
    while BDIO_seek!(fb, 2) == true push!(obs[j], read_uwreal(fb)) end 
    BDIO_close!(fb)
end
mpi = obs[1]
mk = obs[2]
m12 = obs[3] ## not improved, not renormalized
m13 = obs[4]
fpi = obs[5] ## improved & renormalized
fk = obs[6]

phi2 = [8 * t0 * mpi[i] ^ 2 for i in 1:length(mpi)]
phi4 = Array{uwreal,1}()
fpik = Array{uwreal,1}()
c=0
for i in 1:length(mpi)
    for j in 2*i-1+c:2*i+1+c
        push!(phi4, 8 * t0 * (mk[j] ^ 2 + 0.5 * mpi[i] ^ 2))
        push!(fpik, sqrt(t0) * 2/3 * (fk[j] + 0.5 * fpi[i]))
    end
    c+=1
end
uwerr.(phi2)
uwerr.(phi4)
uwerr.(m12)
uwerr.(fpik)

mpi_w, mk_w, fpi_w, fk_w = fve(mpi_w, mk_w, fpi_w, fk_w, ens)
c=0
for i in 1:length(mpi)
    for j in 2*i-1+c:2*i+1+c
        global a, mk[j], b, fk[j] = fve(mpi[i], mk[j], fpi[i], fk[j], ens)
    end
    mpi[i] = a
    fpi[i] = b
    c+=1
end

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
fpi_w_sh = fpi_w + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[24:26]) / sqrt(t0)
fk_w_sh = fk_w + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[27:29]) / sqrt(t0) 
phi2_w_sh = phi2_w + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[4:6])
t0_sh = t0 + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[7:9])
m12_w_sh = m12_w + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[30:34]) * ZP / sqrt(t0) / ZA ## not improved, not renormalized
m12_w_sh_I = (1 + beta_bap[ens.beta] * m12_w_sh) * m12_w_sh
m12_w_I = (1 + beta_bap[ens.beta] * m12_w) * m12_w
m13_w_sh = m13_w + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[35:39]) * ZP / sqrt(t0) / ZA ## not improved, not renormalized
m13_w_sh_I = (1 + beta_bap[ens.beta] * m13_w_sh) * m13_w_sh
m13_w_I = (1 + beta_bap[ens.beta] * m13_w) * m13_w

fpik_sh = [fpik[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[10:12]) for i in 1:length(fpik)]
fpi_sh = [fpi[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[40:42]) for i in 1:length(fpi)]
fk_sh = [fk[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[43:45]) for i in 1:length(fk)]
phi2_sh = [phi2[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[13:15]) for i in 1:length(phi2)]
phi4_sh = [phi4[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[16:18]) for i in 1:length(phi4)]
m12_sh = [m12[i] + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[19:23]) * ZP / sqrt(t0) / ZA for i in 1:length(m12)]
uwerr.(fpik_sh)
uwerr.(fpi_sh)
uwerr.(fk_sh)
uwerr.(phi2_sh)
uwerr.(phi4_sh)
uwerr.(m12_sh)

#======== match & full twist =======#

y = [m12_sh, phi2_sh, phi4_sh]
target = [0 .* m12_sh, [phi2_w_sh for i in 1:length(phi2)], [phi4_ph for i in 1:length(phi4)]]
y = y .- target

kappa, mul, mus = ens_kappa[id], ens_mul[id], ens_mus[id]
kappa_aux = [[kappa[1] for i in 1:9]; [kappa[2] for i in 1:9]; [kappa[3] for i in 1:9]]
mul_aux = [[mul[1] for i in 1:3]; [mul[2] for i in 1:3]; [mul[3] for i in 1:3]]
mus_aux = [mus; mus; mus]
x_l = [[[kappa[1] for i in 1:3]; [kappa[2] for i in 1:3]; [kappa[3] for i in 1:3]] [mul; mul; mul] [mul; mul; mul]]
x_s = [kappa_aux [mul_aux; mul_aux; mul_aux] [mus_aux; mus_aux; mus_aux]]
x = [x_l, x_l, x_s]

#up, chi2, chi_exp, pv = fit_alg([match_m12, match_phi2, match_phi4],x,y,11,wpm=wpm) ##kappa->up[1], mul->up[2], mus->up[3]
if ens.id in ["H102r001", "H102r002"]
    lb = [kappa[1], mul[1], mus[1], -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
    ub = [kappa[2], mul[2], mus[2], Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]
    p0 = [kappa[2], mul[2], mus[2], 0.31773171018397933, -0.000264588017236073, 70.35615274515222, 158.39403084751487, -0.04214173704892017, -0.12252811003970752, 2.712514933784784, 85.2047801509684]
    up, chi2, chi_exp, pv = fit_alg([match_m12, match_phi2, match_phi4],x,y,11,lb,ub,p0,wpm=wpm) ##kappa->up[1], mul->up[2], mus->up[3]
elseif ens.id in ["H105"]
    up, chi2, chi_exp, pv = fit_alg([match_m12, match_phi2, match_phi4],x,y,11,[kappa[2], mul[2], mus[2]],wpm=wpm) ##kappa->up[1], mul->up[2], mus->up[3]
elseif ens.id in ["H105r005"]
    up, chi2, chi_exp, pv = fit_alg_LBFGS([match_m12, match_phi2, match_phi4],x,y,11,[kappa[2], mul[2], mus[2]],wpm=wpm) ##kappa->up[1], mul->up[2], mus->up[3]
else
    up, chi2, chi_exp, pv = fit_alg([match_m12, match_phi2, match_phi4],x,y,11,[kappa[2], mul[2], mus[2]],wpm=wpm) ##kappa->up[1], mul->up[2], mus->up[3]
end

matching_constTr_plot()

#========= interpolate fpik ========#

y = fpik_sh

cc = 0
while cc == 0
    try
    global up_fpik, chi2, chi_exp, pv = fit_alg(interp_fpik_constTr,x_s,y,5,ens_up_fpik[ens.id].+randn(5),wpm=wpm)
    global fpik_matched = interp_fpik_constTr([up[1] up[2] up[3]],up_fpik)[1]
    uwerr(fpik_matched)
    cc+=1
    catch e
    end
end

interp_fpik_constTr_plot()

y = fpi_sh

cc = 0
while cc == 0
    try
    global up_fpi, chi2, chi_exp, pv = fit_alg(interp_fpik_sym,x_l,y,4,ens_up_fpi[ens.id].+randn(4),wpm=wpm)
    global fpi_matched = interp_fpik_sym([up[1] up[2]],up_fpi)[1]
    uwerr(fpi_matched)
    cc+=1
    catch e
    end
end

y = fk_sh

cc = 0
while cc == 0
    try
    global up_fk, chi2, chi_exp, pv = fit_alg(interp_fpik_constTr,x_s,y,5,ens_up_fk[ens.id].+randn(5),wpm=wpm)
    global fk_matched = interp_fpik_constTr([up[1] up[2] up[3]],up_fk)[1]
    uwerr(fk_matched)
    cc+=1
    catch e
    end
end


#========= save bdio ===============#

obs = [t0_sh, phi2_w_sh, m12_w_sh_I, m13_w_sh_I, fpi_w_sh, fk_w_sh, fpik_w_sh]
fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted/", ens.id, "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

obs = [up[1], up[2], up[3], fpi_matched, fk_matched, fpik_matched]
fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted/", ens.id, "_obs_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)









