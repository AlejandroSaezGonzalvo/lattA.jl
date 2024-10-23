#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, PyCall, LsqFit, LinearAlgebra
using ADerrors: err
using lattA: fit_alg

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/chiral-continuum_fits.jl");

new = true
if new == true
    ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "D450", "N202", "N203", "D200", "E250", "N300", "N302", "J303", "E300", "J500", "J501"]
    ens_old = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "D200", "N300", "J303", "J501"]
    ens_new = ["D450", "E250", "N302", "E300", "J500", "J501"]
    ens_av = ["H101", "H102", "H105", "H400", "D450", "N202", "N203", "D200", "E250", "N300", "N302", "J303", "E300", "J500", "J501"]
else
    ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "D200", "N300", "J303"]
    ens_old = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "D200", "N300", "J303"]
    ens_new = [""]
    ens_av = ["H101", "H102", "H105", "H400", "N202", "N203", "D200", "N300", "J303"]
end
ens_sym = ["H101", "H400", "N202", "N300", "J500"]
ens_nosym = ["H102", "H105", "D450", "N203", "D200", "E250", "N302", "J303", "E300", "J501"]
ens_nosym355 = ["N203", "D200", "E250", "N302", "J303", "E300", "J501"]
ens_40 = ["H101", "H102", "H400", "D450", "N202", "N203", "D200", "E250", "N300", "N302", "J303", "E300", "J500", "J501"]
ens_41 = ["H101", "H102", "H400", "D450", "N202", "N203", "D200", "N300", "N302", "J303", "E300", "J500", "J501"]
ens_42 = ["H101", "H102", "H400", "D450", "N202", "N203", "D200", "N300", "N302", "E300", "J500", "J501"]
ens_40_nosym = ["H102", "D450", "N203", "D200", "E250", "N302", "J303", "E300", "J501"]
ens_41_nosym = ["H102", "D450", "N203", "D200", "N302", "J303", "E300", "J501"]
ens_42_nosym = ["H102", "D450", "N203", "D200", "N302", "E300", "J501"]
ind_sym = findall(x -> x in ens_sym, ens_av)
ind_nosym = findall(x -> x in ens_nosym, ens_av)
ind_nosym355 = findall(x -> x in ens_nosym355, ens_av)
ens_340 = findall(x -> x in ["H101", "H102", "H105"], ens_av)
ens_346 = findall(x -> x in ["H400", "D450"], ens_av)
ens_355 = findall(x -> x in ["N202", "N203", "D200", "E250"], ens_av)
ens_370 = findall(x -> x in ["N300", "N302", "J303", "E300"], ens_av)
ens_385 = findall(x -> x in ["J500", "J501"], ens_av)
ind_phi204 = findall(x -> x in ["H105", "H105r005", "D450", "D200", "E250", "J303", "E300"], ens_av)
ind_mL_40 = findall(x -> x in ens_40, ens_av)
ind_mL_41 = findall(x -> x in ens_41, ens_av)
ind_mL_42 = findall(x -> x in ens_42, ens_av)
ind_phi205 = [2,7,12,16]
ens_340_nosym = findall(x -> x in ["H102", "H105"], ens_av)
ens_346_nosym = findall(x -> x in ["D450"], ens_av)
ens_355_nosym = findall(x -> x in ["N203", "D200", "E250"], ens_av)
ens_370_nosym = findall(x -> x in ["N302", "J303", "E300"], ens_av)
ens_385_nosym = findall(x -> x in ["J501"], ens_av)
ind_mL_40_nosym = findall(x -> x in ens_40_nosym, ens_av)
ind_mL_41_nosym = findall(x -> x in ens_41_nosym, ens_av)
ind_mL_42_nosym = findall(x -> x in ens_42_nosym, ens_av)


beta = [3.4,3.4,3.4,3.46,3.46,3.55,3.55,3.55,3.55,3.55,3.7,3.7,3.7,3.7,3.85,3.85]

#============================== read BDIO =====================================================================================#

    obs = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_tm_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    m13_tm_sh = Array{uwreal,1}()
    for i in 1:length(ens)
        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted_FLAG21/", ens[i], "_m13_matched_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
        BDIO_seek!(fb); push!(m13_tm_sh, read_uwreal(fb))
        BDIO_close!(fb)

        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ens[i], "_obs_wil_un.bdio"), "r")
        BDIO_seek!(fb); push!(obs[i], read_uwreal(fb))
        while BDIO_seek!(fb, 2) == true push!(obs[i], read_uwreal(fb)) end
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
    t0 = [obs[i][1] for i in 1:length(obs)]
    mpi = [obs[i][2] for i in 1:length(obs)]
    mk = [obs[i][3] for i in 1:length(obs)]
    m12 = [obs[i][4] for i in 1:length(obs)]
    m13 = [obs[i][5] for i in 1:length(obs)]
    fpi = [obs[i][6] for i in 1:length(obs)]
    fk = [obs[i][7] for i in 1:length(obs)]
    t0_sh = [obs_sh[i][1] for i in 1:length(obs)]
    phi2_sh = [obs_sh[i][2] for i in 1:length(obs)]
    m12_sh = [obs_sh[i][3] for i in 1:length(obs)]
    m13_sh = [obs_sh[i][4] for i in 1:length(obs)]
    fpi_sh = [obs_sh[i][5] for i in 1:length(obs)]
    fk_sh = [obs_sh[i][6] for i in 1:length(obs)]
    t0fpik_st_sh = [obs_sh[i][7] for i in 1:length(obs)]
    kappa = [obs_tm_sh[i][1] for i in 1:length(obs_tm_sh)]
    mul = [obs_tm_sh[i][2] for i in 1:length(obs_tm_sh)]
    mus = [obs_tm_sh[i][3] for i in 1:length(obs_tm_sh)]
    fpi_matched = [obs_tm_sh[i][4] for i in 1:length(obs_tm_sh)]
    fk_matched = [obs_tm_sh[i][5] for i in 1:length(obs_tm_sh)]
    t0fpik_sh = [obs_tm_sh[i][6] for i in 1:length(obs_tm_sh)]

    ind = ensemble_inv["H102r002"]
    t0[ind-1] = plat_av(t0, [ind-1,ind]); deleteat!(t0, ind)
    mpi[ind-1] = plat_av(mpi, [ind-1,ind]); deleteat!(mpi, ind)
    mk[ind-1] = plat_av(mk, [ind-1,ind]); deleteat!(mk, ind)
    m12[ind-1] = plat_av(m12, [ind-1,ind]); deleteat!(m12, ind)
    m13[ind-1] = plat_av(m13, [ind-1,ind]); deleteat!(m13, ind)
    fpi[ind-1] = plat_av(fpi, [ind-1,ind]); deleteat!(fpi, ind)
    fk[ind-1] = plat_av(fk, [ind-1,ind]); deleteat!(fk, ind)
    t0_sh[ind-1] = plat_av(t0_sh, [ind-1,ind]); deleteat!(t0_sh, ind)
    phi2_sh[ind-1] = plat_av(phi2_sh, [ind-1,ind]); deleteat!(phi2_sh, ind)
    m12_sh[ind-1] = plat_av(m12_sh, [ind-1,ind]); deleteat!(m12_sh, ind)
    m13_sh[ind-1] = plat_av(m13_sh, [ind-1,ind]); deleteat!(m13_sh, ind)
    mul[ind-1] = plat_av(mul, [ind-1,ind]); deleteat!(mul, ind)
    mus[ind-1] = plat_av(mus, [ind-1,ind]); deleteat!(mus, ind)
    fpi_sh[ind-1] = plat_av(fpi_sh, [ind-1,ind]); deleteat!(fpi_sh, ind)
    fk_sh[ind-1] = plat_av(fk_sh, [ind-1,ind]); deleteat!(fk_sh, ind)
    t0fpik_st_sh[ind-1] = plat_av(t0fpik_st_sh, [ind-1,ind]); deleteat!(t0fpik_st_sh, ind)
    fpi_matched[ind-1] = plat_av(fpi_matched, [ind-1,ind]); deleteat!(fpi_matched, ind)
    fk_matched[ind-1] = plat_av(fk_matched, [ind-1,ind]); deleteat!(fk_matched, ind)
    t0fpik_sh[ind-1] = plat_av(t0fpik_sh, [ind-1,ind]); deleteat!(t0fpik_sh, ind)
    m13_tm_sh[ind-1] = plat_av(m13_tm_sh, [ind-1,ind]); deleteat!(m13_tm_sh, ind)

    ind = ensemble_inv["H105r005"] - 1
    t0[ind-1] = plat_av(t0, [ind-1,ind]); deleteat!(t0, ind)
    mpi[ind-1] = plat_av(mpi, [ind-1,ind]); deleteat!(mpi, ind)
    mk[ind-1] = plat_av(mk, [ind-1,ind]); deleteat!(mk, ind)
    m12[ind-1] = plat_av(m12, [ind-1,ind]); deleteat!(m12, ind)
    m13[ind-1] = plat_av(m13, [ind-1,ind]); deleteat!(m13, ind)
    fpi[ind-1] = plat_av(fpi, [ind-1,ind]); deleteat!(fpi, ind)
    fk[ind-1] = plat_av(fk, [ind-1,ind]); deleteat!(fk, ind)
    t0_sh[ind-1] = plat_av(t0_sh, [ind-1,ind]); deleteat!(t0_sh, ind)
    phi2_sh[ind-1] = plat_av(phi2_sh, [ind-1,ind]); deleteat!(phi2_sh, ind)
    m12_sh[ind-1] = plat_av(m12_sh, [ind-1,ind]); deleteat!(m12_sh, ind)
    m13_sh[ind-1] = plat_av(m13_sh, [ind-1,ind]); deleteat!(m13_sh, ind)
    mul[ind-1] = plat_av(mul, [ind-1,ind]); deleteat!(mul, ind)
    mus[ind-1] = plat_av(mus, [ind-1,ind]); deleteat!(mus, ind)
    fpi_sh[ind-1] = plat_av(fpi_sh, [ind-1,ind]); deleteat!(fpi_sh, ind)
    fk_sh[ind-1] = plat_av(fk_sh, [ind-1,ind]); deleteat!(fk_sh, ind)
    t0fpik_st_sh[ind-1] = plat_av(t0fpik_st_sh, [ind-1,ind]); deleteat!(t0fpik_st_sh, ind)
    fpi_matched[ind-1] = plat_av(fpi_matched, [ind-1,ind]); deleteat!(fpi_matched, ind)
    fk_matched[ind-1] = plat_av(fk_matched, [ind-1,ind]); deleteat!(fk_matched, ind)
    t0fpik_sh[ind-1] = plat_av(t0fpik_sh, [ind-1,ind]); deleteat!(t0fpik_sh, ind)
    m13_tm_sh[ind-1] = plat_av(m13_tm_sh, [ind-1,ind]); deleteat!(m13_tm_sh, ind)

    uwerr.(phi2_sh)
    phi4_sh = [phi4_ph for i in 1:length(phi2_sh)]
    phi2_sym = 2 / 3 * phi4_sh
    phi2_sym_ph = 2 / 3 * phi4_ph
    t0fpi_st_sh = sqrt.(8 * t0_sh) .* fpi_sh
    t0fk_st_sh = sqrt.(8 * t0_sh) .* fk_sh
    t0fpi_sh = sqrt.(8 * t0_sh) .* fpi_matched
    t0fk_sh = sqrt.(8 * t0_sh) .* fk_matched
    t0fpik_sh = sqrt(8) .* t0fpik_sh
    t0fpik_st_sh = sqrt(8) .* t0fpik_st_sh

    ZA = [beta_ZA[beta[i]] for i in 1:length(ens_av)]
    rat = ZA .* 2 .* m13_tm_sh ./ mus; uwerr.(rat)

#==============================================================================================================================#

#========================== plot =============================#

    errorbar(value.(phi2_sh[ens_340]), value.(rat[ens_340]), err.(rat[ens_340]), fmt="x", color="purple")
    errorbar(value.(phi2_sh[ens_346]), value.(rat[ens_346]), err.(rat[ens_346]), fmt="x", color="green")
    errorbar(value.(phi2_sh[ens_355]), value.(rat[ens_355]), err.(rat[ens_355]), fmt="x", color="blue")
    errorbar(value.(phi2_sh[ens_370]), value.(rat[ens_370]), err.(rat[ens_370]), fmt="x", color="orange")
    errorbar(value.(phi2_sh[ens_385]), value.(rat[ens_385]), err.(rat[ens_385]), fmt="x", color="red")
    ylabel(L"$2Z_Am_{13}^{Wtm}/\mu_s$")
    xlabel(L"$\phi_2$")
    tight_layout()

#=============================================================#

#========================= fit ===============================#

    function f(x,p)
        return [p[1] + p[2] * x[i,2] + p[3] * x[i,2] ^ 2 + p[4] * x[i,1] for i in 1:length(x[:,1])]
    end

    x = [1 ./ value.(t0_sh[ind_nosym]) value.(phi2_sh[ind_nosym])]
    y = rat[ind_nosym]
    uwerr.(y)

    up, chi2, chi_exp, pv = fit_alg(f,x,y,4)

    errorbar(value.(phi2_sh[ens_340]), value.(rat[ens_340]), err.(rat[ens_340]), fmt="x", color="purple")
    errorbar(value.(phi2_sh[ens_346]), value.(rat[ens_346]), err.(rat[ens_346]), fmt="x", color="green")
    errorbar(value.(phi2_sh[ens_355]), value.(rat[ens_355]), err.(rat[ens_355]), fmt="x", color="blue")
    errorbar(value.(phi2_sh[ens_370]), value.(rat[ens_370]), err.(rat[ens_370]), fmt="x", color="orange")
    errorbar(value.(phi2_sh[ens_385]), value.(rat[ens_385]), err.(rat[ens_385]), fmt="x", color="red")

    

    ylabel(L"$2Z_Am_{13}^{Wtm}/\mu_s$")
    xlabel(L"$\phi_2$")
    tight_layout()

#=============================================================#