#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, LsqFit, LinearAlgebra
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/chiral-continuum_fits.jl");

ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303", "E300", "J500"]
ens_old = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "N200", "D200", "N300", "J303"]
ens_new = ["D450", "E250", "N302", "E300", "J500"]
ens_av = ["H101", "H102", "H105", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303", "E300", "J500"]
ens_sym = ["H101", "H400", "N202", "N300", "J500"]
ens_nosym = ["H102", "H105", "D450", "N203", "N200", "D200", "E250", "N302", "J303", "E300"]
ens_nosym355 = ["N203", "N200", "D200", "E250", "N302", "J303", "E300"]
ens_40 = ["H101", "H102", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303", "E300", "J500"]
ens_41 = ["H101", "H102", "H400", "D450", "N202", "N203", "N200", "D200", "N300", "N302", "J303", "E300", "J500"]
ens_42 = ["H101", "H102", "H400", "D450", "N202", "N203", "N200", "D200", "N300", "N302", "E300", "J500"]
ens_340 = findall(x -> x in ["H101", "H102", "H105"], ens_av)
ens_346 = findall(x -> x in ["H400", "D450"], ens_av)
ens_355 = findall(x -> x in ["N202", "N203", "N200", "D200", "E250"], ens_av)
ens_370 = findall(x -> x in ["N300", "N302", "J303", "E300"], ens_av)
ens_385 = findall(x -> x in ["J500"], ens_av)
ind_sym = findall(x -> x in ens_sym, ens_av)
ind_nosym = findall(x -> x in ens_nosym, ens_av)
ind_nosym355 = findall(x -> x in ens_nosym355, ens_av)
ind_phi204 = findall(x -> x in ["H105", "H105r005", "D450", "N200", "D200", "E250", "J303", "E300"], ens_av)
ind_mL_40 = findall(x -> x in ens_40, ens_av)
ind_mL_41 = findall(x -> x in ens_41, ens_av)
ind_mL_42 = findall(x -> x in ens_42, ens_av)

#============================== read BDIO =====================================================================================#

    obs = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_tm_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    for i in 1:length(ens)
        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ens[i], "_obs_wil_un.bdio"), "r")
        BDIO_seek!(fb); push!(obs[i], read_uwreal(fb))
        while BDIO_seek!(fb, 2) == true push!(obs[i], read_uwreal(fb)) end
        BDIO_close!(fb)

        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted/", ens[i], "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
        BDIO_seek!(fb); push!(obs_sh[i], read_uwreal(fb))
        while BDIO_seek!(fb, 2) == true push!(obs_sh[i], read_uwreal(fb)) end
        BDIO_close!(fb)

        fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/shifted/", ens[i], "_obs_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
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
    m12_sh = [obs_sh[i][3] * beta_ZA[EnsInfo(ens[i], ens_db[ens[i]]).beta] / beta_ZP[EnsInfo(ens[i], ens_db[ens[i]]).beta] for i in 1:length(obs)]
    m13_sh = [obs_sh[i][4] * beta_ZA[EnsInfo(ens[i], ens_db[ens[i]]).beta] / beta_ZP[EnsInfo(ens[i], ens_db[ens[i]]).beta] for i in 1:length(obs)]
    fpi_sh = [obs_sh[i][5] for i in 1:length(obs)]
    fk_sh = [obs_sh[i][6] for i in 1:length(obs)]
    t0fpik_st_sh = [obs_sh[i][7] for i in 1:length(obs)]
    kappa = [obs_tm_sh[i][1] for i in 1:length(obs_tm_sh)]
    mul = [obs_tm_sh[i][2] / beta_ZP[EnsInfo(ens[i], ens_db[ens[i]]).beta] for i in 1:length(obs_tm_sh)]
    mus = [obs_tm_sh[i][3] / beta_ZP[EnsInfo(ens[i], ens_db[ens[i]]).beta] for i in 1:length(obs_tm_sh)]
    fpi_matched = [obs_tm_sh[i][4] for i in 1:length(obs_tm_sh)]
    fk_matched = [obs_tm_sh[i][5] for i in 1:length(obs_tm_sh)]
    t0fpik_sh = [obs_tm_sh[i][6] for i in 1:length(obs_tm_sh)]

    ind = ensemble_inv["H102r002"]
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

    ind = ensemble_inv["H105r005"] - 1
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
    phi12_st = sqrt.(8 * t0_sh) .* m12_sh
    phi13_st = sqrt.(8 * t0_sh) .* m13_sh
    phi12 = sqrt.(8 * t0_sh) .* mul
    phi13 = sqrt.(8 * t0_sh) .* (mul .+ mus) ./ 2
    phiK_sh = phi4_sh .- 0.5 * phi2_sh

    y1_st = phi12_st ./ phi13_st; uwerr.(y1_st)
    y1 = phi12 ./ phi13; uwerr.(y1)
    y2_st = 2 * phi13_st ./ phiK_sh + phi12_st ./ phi2_sh; uwerr.(y2_st)
    y2 = 2 * phi13 ./ phiK_sh + phi12 ./ phi2_sh; uwerr.(y2)

    #t0fpik_sh = sqrt.(8 * t0_sh) .* (2/3) .* (1/2 * fpi_matched .+ fk_matched)
    #t0fpik_st_sh = sqrt.(8 * t0_sh) .* (2/3) .* (1/2 * fpi_sh .+ fk_sh) 

    t0_sh_sym = [[t0_sh[1] for i in 1:3]; [t0_sh[4]]; [t0_sh[5] for i in 5:8]; [t0_sh[9] for i in 9:10]]
    #t0fpik_sh = t0fpik_sh ./ sqrt.(t0_sh) .* sqrt.(t0_sh_sym)
    #t0fpik_st_sh = t0fpik_st_sh ./ sqrt.(t0_sh) .* sqrt.(t0_sh_sym)

#==============================================================================================================================#

#============================== plots =========================================================================================#

    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\phi_{12}/\phi_{13}$")
    errorbar(value.(phi2_sh[ens_340]), value.(y1[ens_340]), err.(y1[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(y1[ens_346]), err.(y1[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(y1[ens_355]), err.(y1[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(y1[ens_370]), err.(y1[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(y1[ens_385]), err.(y1[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    #ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/phi12_div_phi13_tm.pdf")

    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\phi_{12}/\phi_{13}$")
    errorbar(value.(phi2_sh[ens_340]), value.(y1_st[ens_340]), err.(y1_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(y1_st[ens_346]), err.(y1_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(y1_st[ens_355]), err.(y1_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(y1_st[ens_370]), err.(y1_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(y1_st[ens_385]), err.(y1_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    #ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/phi12_div_phi13_w.pdf")

    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$2\phi_{13}/\phi_K+\phi_{12}/\phi_2$")
    errorbar(value.(phi2_sh[ens_340]), value.(y2[ens_340]), err.(y2[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(y2[ens_346]), err.(y2[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(y2[ens_355]), err.(y2[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(y2[ens_370]), err.(y2[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(y2[ens_385]), err.(y2[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    #ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/2phi13_div_phiK_plus_phi12_div_phi2_tm.pdf")

    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$2\phi_{13}/\phi_K+\phi_{12}/\phi_2$")
    errorbar(value.(phi2_sh[ens_340]), value.(y2_st[ens_340]), err.(y2_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(y2_st[ens_346]), err.(y2_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(y2_st[ens_355]), err.(y2_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(y2_st[ens_370]), err.(y2_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(y2_st[ens_385]), err.(y2_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    #ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/2phi13_div_phiK_plus_phi12_div_phi2_w.pdf")

#==============================================================================================================================#

#============================== chiral & continuum limits =====================================================================#

    ## prepare data SU3 & Tay
        TIC = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        chi = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        W = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        y1_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        y2_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        pval = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        uprm_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]

        x = [1 ./ (8 .* t0_sh) phi2_sh phi4_sh phi2_sym]
        x_combined = [x; x]
        x_355 = x[[ens_346; ens_355; ens_370; ens_385],:]
        x_combined_355 = [x_355; x_355]
        x_combined_355_only_Wil = [x_355; x]
        x_combined_355_only_Wtm = [x; x_355]
        x_nosym = x[ind_nosym,:]
        x_nosym355 = x[ind_nosym355,:]
        x_combined_nosym = [x_nosym; x_nosym]
        x_combined_nosym355 = [x_nosym355; x_nosym355]
        x_mL_40 = x[ind_mL_40, :]
        x_mL_41 = x[ind_mL_41, :]
        x_mL_42 = x[ind_mL_42, :]
        x_combined_mL_40 = [x_mL_40; x_mL_40]
        x_combined_mL_41 = [x_mL_41; x_mL_41]
        x_combined_mL_42 = [x_mL_42; x_mL_42]
        x_ph = [0.0 (phi2_ph) (phi4_ph) (phi2_sym_ph)]
        x_3552 = x[[ens_355; ens_370; ens_385],:]
        x_combined_3552 = [x_3552; x_3552]
        x_phi204 = x[ind_phi204,:]
        x_combined_phi204 = [x_phi204; x_phi204]

        y = [y1; y2]
        y_st = [y1_st; y2_st]
        y_combined = [y_st; y]
        y_355 = [y1[[ens_346; ens_355; ens_370; ens_385]]; y2[[ens_346; ens_355; ens_370; ens_385]]]
        y_st_355 = [y1_st[[ens_346; ens_355; ens_370; ens_385]]; y2_st[[ens_346; ens_355; ens_370; ens_385]]]
        y_combined_355 = [y_st_355; y_355]
        y_3552 = [y1[[ens_355; ens_370; ens_385]]; y2[[ens_355; ens_370; ens_385]]]
        y_st_3552 = [y1_st[[ens_355; ens_370; ens_385]]; y2_st[[ens_355; ens_370; ens_385]]]
        y_combined_3552 = [y_st_3552; y_3552]
        y_nosym = [y1[ind_nosym]; y2[ind_nosym]]
        y_st_nosym = [y1_st[ind_nosym]; y2_st[ind_nosym]]
        y_combined_nosym = [y_st_nosym; y_nosym]
        y_nosym355 = [y1[ind_nosym355]; y2[ind_nosym355]]
        y_st_nosym355 = [y_st[ind_nosym355]; y2_st[ind_nosym355]]
        y_combined_nosym355 = [y_st_nosym355; y_nosym355]
        y_phi204 = [y1[ind_phi204]; y2[ind_phi204]]
        y_st_phi204 = [y1_st[ind_phi204]; y2_st[ind_phi204]]
        y_combined_phi204 = [y_st_phi204; y_phi204]
        y_mL_42 = [y1[ind_mL_42]; y2[ind_mL_42]]
        y_st_mL_42 = [y1_st[ind_mL_42]; y2_st[ind_mL_42]]
        y_combined_mL_42 = [y_st_mL_42; y_mL_42]

        cuts_y = [y, y_355, y_3552, y_nosym, y_nosym355, y_phi204, y_mL_42]
        cuts_y_st = [y_st, y_st_355, y_st_3552, y_st_nosym, y_st_nosym355, y_st_phi204, y_st_mL_42]
        cuts_y_combined = [y_combined, y_combined_355, y_combined_3552, y_combined_nosym, y_combined_nosym355, y_combined_phi204, y_combined_mL_42]
        cuts_x = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_st = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_combined = [x_combined, x_combined_355, x_combined_3552, x_combined_nosym, x_combined_nosym355, x_combined_phi204, x_combined_mL_42]

        set_y = [cuts_y, cuts_y_st, cuts_y_combined]
        set_x = [cuts_x, cuts_x_st, cuts_x_combined]

        Wm = [inv.(Symmetric.(cov.(set_y[i]))) for i in 1:length(set_y)]
        Wm = [convert.(Matrix{Float64}, Wm[i]) for i in 1:length(set_y)]
    ##

    models = [mR_a2; mR_aas; mR_a2phi2; mR_Tay_a2; mR_Tay_aas; mR_Tay_a2phi2]
    models_combined = [mR_a2_combined; mR_aas_combined; mR_a2phi2_combined; mR_a2a2phi2_combined; mR_a2phi2a2_combined; mR_Tay_a2_combined; mR_Tay_aas_combined; mR_Tay_a2phi2_combined; mR_Tay_a2a2phi2_combined; mR_Tay_a2phi2a2_combined]
    models = [models, models, models_combined]
    param = [6,6,8,7,7,9]
    param_combined = [8,8,12,10,10,9,9,13,11,11]
    param = [param, param, param_combined]

    for k in 1:length(set_y)
        for i in 1:length(models[k])
            for j in 1:length(cuts_y)
                x = set_x[k][j]
                y = set_y[k][j]
                global L1 = length(set_y[1][j])
                global L2 = length(set_y[1][j])
                uprm, chi_exp, chi2, pval_aux, doff = fit_alg(models[k][i], value.([x; x]), y, param[k][i])
                push!(TIC[k], chi2 - 2 * chi_exp)
                push!(pval[k], pval_aux)
                push!(y1_ph_vec[k], models[k][i]([x_ph;x],uprm)[1])
                push!(y2_ph_vec[k], models[k][i]([x_ph;x],uprm)[2])
                if i == j == 1 uprm_plot[k] = uprm end
            end
        end
    end

    TIC = [TIC[k] .- minimum.(TIC)[k] for k in 1:length(TIC)]
    W = [exp.(-0.5 * TIC[k]) ./ sum(exp.(-0.5 * TIC[k])) for k in 1:length(TIC)]
    y1_ph = [sum(y1_ph_vec[k] .* W[k]) for k in 1:length(TIC)]
    syst1 = [sqrt(sum(y1_ph_vec[k] .^ 2 .* W[k]) - (sum(y1_ph_vec[k] .* W[k])) ^ 2) for k in 1:length(TIC)]
    y1_ph = y1_ph .+ [uwreal([0.0, value(syst1[k])], string("syst mR chiral ", k, " 3rd")) for k in 1:length(TIC)]
    uwerr.(y1_ph)
    y2_ph = [sum(y2_ph_vec[k] .* W[k]) for k in 1:length(TIC)]
    syst2 = [sqrt(sum(y2_ph_vec[k] .^ 2 .* W[k]) - (sum(y2_ph_vec[k] .* W[k])) ^ 2) for k in 1:length(TIC)]
    y2_ph = y2_ph .+ [uwreal([0.0, value(syst2[k])], string("syst mR chiral ", k, " 3rd")) for k in 1:length(TIC)]
    uwerr.(y2_ph)

    W_aux = deepcopy(W)
    pval_aux = deepcopy(pval)
    [uwerr.(y1_ph_vec[k]) for k in 1:length(TIC)]
    [uwerr.(y2_ph_vec[k]) for k in 1:length(TIC)]

    
