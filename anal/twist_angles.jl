#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot
using ADerrors: err
using lattA: EnsInfo, fve, fit_alg

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");

ens = ["H102r001", "H102r002", "N203", "N302", "J501"]
ens_av = ["H102", "N203", "N302", "J501"]
beta = [3.4, 3.55, 3.7, 3.85]

#================== read BDIO ===================================#

    m13 = Array{uwreal,1}()
    obs_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_tm_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    for i in 1:length(ens)
        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted_FLAG21/", ens[i], "_m13_matched_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
        BDIO_seek!(fb); push!(m13, read_uwreal(fb))
        BDIO_close!(fb)

        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted_FLAG21/", ens[i], "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
        BDIO_seek!(fb); push!(obs_sh[i], read_uwreal(fb))
        while BDIO_seek!(fb, 2) == true push!(obs_sh[i], read_uwreal(fb)) end
        BDIO_close!(fb)

        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted_FLAG21/", ens[i], "_obs_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
        BDIO_seek!(fb); push!(obs_tm_sh[i], read_uwreal(fb))
        while BDIO_seek!(fb, 2) == true push!(obs_tm_sh[i], read_uwreal(fb)) end
        BDIO_close!(fb)
    end
    t0 = [obs_sh[i][1] for i in 1:length(obs_sh)]
    mus = [obs_tm_sh[i][3] for i in 1:length(obs_tm_sh)]
    phi2 = [obs_sh[i][2] for i in 1:length(obs_sh)]    

    t0[1] = plat_av(t0, [1,2]); deleteat!(t0, 2)
    mus[1] = plat_av(mus, [1,2]); deleteat!(mus, 2)
    m13[1] = plat_av(m13, [1,2]); deleteat!(m13, 2)
    phi2[1] = plat_av(phi2, [1,2]); deleteat!(phi2, 2)

    uwerr.(mus)
    uwerr.(m13)
    uwerr.(phi2)
    phi4_sh = [phi4_ph for i in 1:length(phi2_sh)]

    ZA = [beta_ZA[beta[i]] for i in 1:length(ens_av)]
    ZP = [beta_ZP[beta[i]] for i in 1:length(ens_av)]
    m13_R = ZA ./ ZP .* m13
    mus_R = 1 ./ ZP .* mus
    rat = m13_R ./ mus_R; uwerr.(rat)

#================================================================#

#================== plot ========================================#

    errorbar(value.(1 ./ t0), value.(rat), err.(rat), fmt="x", label=L"$\phi_2\approx0.5, \phi_4=1.096(11)$")
    xlabel(L"$a^2/t_0$")
    ylabel(L"$am_{13}^{(v), R}/a\mu_s^{(v), R}$")
    legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/twist_angle_a2.pdf")

#================================================================#


