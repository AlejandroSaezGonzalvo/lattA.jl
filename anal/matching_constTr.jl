import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H105"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#=========== read bdio ================#

obs = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "r")
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
    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "r")
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
phi2_w_sh = phi2_w + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[4:6])
t0_sh = t0 + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[7:9])
#m12_w_sh = m12_w + (phi4_ph - phi4_w) * md_?(phi2_w, 1/t0, par[???]) * ZP / sqrt(t0) ## not improved, not renormalized
#m12_w_sh = (1 + beta_bap[ens.beta] * m12_w_sh) * m12_w_sh
#m12_w = (1 + beta_bap[ens.beta] * m12_w) * m12_w

fpik_sh = [fpik[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[10:12]) for i in 1:length(fpik)]
phi2_sh = [phi2[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[13:15]) for i in 1:length(phi2)]
phi4_sh = [phi4[i] + (phi4_ph - phi4_w) * md_1(phi2_w, 1/t0, par[16:18]) for i in 1:length(phi4)]
m12_sh = [m12[i] + (phi4_ph - phi4_w) * md_2(phi2_w, 1/t0, par[19:23]) * ZP / sqrt(t0) for i in 1:length(m12)]
uwerr.(fpik_sh)
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

up, chi2, chi_exp, pv = fit_alg([match_m12, match_phi2, match_phi4],x,y,11,[kappa[2], mul[2], mus[2]],wpm=wpm) ##kappa->up[1], mul->up[2], mus->up[3]

matching_constTr_plot()

#========= interpolate fpik ========#

y = fpik_sh
up_fpik, chi2, chi_exp, pv = fit_alg(interp_fpik_constTr,x_s,y,5,ens_up_fpik[ens.id],wpm=wpm) 
fpik_matched = interp_fpik_constTr([up[1] up[2] up[3]],up_fpik)[1]
uwerr(fpik_matched)

interp_fpik_constTr_plot()

#========= save bdio ===============#

obs = [t0_sh, phi2_w_sh, fpik_w_sh]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

obs = [up[1], up[2], up[3], fpik_matched]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)









