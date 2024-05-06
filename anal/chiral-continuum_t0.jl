#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, LsqFit, LinearAlgebra
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/chiral-continuum_fits.jl");

new = true
if new == true
    ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303", "E300", "J500", "J501"]
    ens_old = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "N200", "D200", "N300", "J303", "J501"]
    ens_new = ["D450", "E250", "N302", "E300", "J500", "J501"]
    ens_av = ["H101", "H102", "H105", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303", "E300", "J500", "J501"]
else
    ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "N200", "D200", "N300", "J303"]
    ens_old = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "N202", "N203", "N200", "D200", "N300", "J303"]
    ens_new = [""]
    ens_av = ["H101", "H102", "H105", "H400", "N202", "N203", "N200", "D200", "N300", "J303"]
end
ens_sym = ["H101", "H400", "N202", "N300", "J500"]
ens_nosym = ["H102", "H105", "D450", "N203", "N200", "D200", "E250", "N302", "J303", "E300", "J501"]
ens_nosym355 = ["N203", "N200", "D200", "E250", "N302", "J303", "E300", "J501"]
ens_40 = ["H101", "H102", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303", "E300", "J500", "J501"]
ens_41 = ["H101", "H102", "H400", "D450", "N202", "N203", "N200", "D200", "N300", "N302", "J303", "E300", "J500", "J501"]
ens_42 = ["H101", "H102", "H400", "D450", "N202", "N203", "N200", "D200", "N300", "N302", "E300", "J500", "J501"]
ens_340 = findall(x -> x in ["H101", "H102", "H105"], ens_av)
ens_346 = findall(x -> x in ["H400", "D450"], ens_av)
ens_355 = findall(x -> x in ["N202", "N203", "N200", "D200", "E250"], ens_av)
ens_370 = findall(x -> x in ["N300", "N302", "J303", "E300"], ens_av)
ens_385 = findall(x -> x in ["J500", "J501"], ens_av)#, "J501"]
ind_sym = findall(x -> x in ens_sym, ens_av)
ind_nosym = findall(x -> x in ens_nosym, ens_av)
ind_nosym355 = findall(x -> x in ens_nosym355, ens_av)
ind_phi204 = findall(x -> x in ["H105", "H105r005", "D450", "N200", "D200", "E250", "J303", "E300"], ens_av)
ind_mL_40 = findall(x -> x in ens_40, ens_av)
ind_mL_41 = findall(x -> x in ens_41, ens_av)
ind_mL_42 = findall(x -> x in ens_42, ens_av)

fpik_add = true

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

    if fpik_add == true
        t0fpik_sh = sqrt.(8 * t0_sh) .* (2/3) .* (1/2 * fpi_matched .+ fk_matched)
        t0fpik_st_sh = sqrt.(8 * t0_sh) .* (2/3) .* (1/2 * fpi_sh .+ fk_sh) 
    end

    t0_sh_sym = [[t0_sh[1] for i in 1:3]; [t0_sh[4] for i in 4:5]; [t0_sh[6] for i in 6:10]; [t0_sh[11] for i in 11:14]; [t0_sh[15] for i in 15:16]]
    #t0fpik_sh = t0fpik_sh ./ sqrt.(t0_sh) .* sqrt.(t0_sh_sym)
    #t0fpik_st_sh = t0fpik_st_sh ./ sqrt.(t0_sh) .* sqrt.(t0_sh_sym)

#==============================================================================================================================#

#============================== plots =========================================================================================#

    uwerr.(t0fpik_sh)
    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{8t_0}f_{\pi K}$")
    errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_sh[ens_340]), err.(t0fpik_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_sh[ens_346]), err.(t0fpik_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_sh[ens_355]), err.(t0fpik_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_sh[ens_370]), err.(t0fpik_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_sh[ens_385]), err.(t0fpik_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0fpik_tm.pdf")

    uwerr.(t0fpik_st_sh)
    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{8t_0}f_{\pi K}$")
    errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_st_sh[ens_340]), err.(t0fpik_st_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_st_sh[ens_346]), err.(t0fpik_st_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_st_sh[ens_355]), err.(t0fpik_st_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_st_sh[ens_370]), err.(t0fpik_st_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_st_sh[ens_385]), err.(t0fpik_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0fpik_st.pdf")

    #t0fpi_st_sh = sqrt.(t0_sh) .* fpi_sh
    uwerr.(t0fpi_st_sh)
    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{8t_0}f_{\pi}$")
    errorbar(value.(phi2_sh[ens_340]), value.(t0fpi_st_sh[ens_340]), err.(t0fpi_st_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(t0fpi_st_sh[ens_346]), err.(t0fpi_st_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(t0fpi_st_sh[ens_355]), err.(t0fpi_st_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(t0fpi_st_sh[ens_370]), err.(t0fpi_st_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(t0fpi_st_sh[ens_385]), err.(t0fpi_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    ax[:set_ylim]([0.24, 0.321])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0fpi_st.pdf")

    #t0fk_st_sh = sqrt.(t0_sh) .* fk_sh
    uwerr.(t0fk_st_sh)
    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{8t_0}f_{K}$")
    errorbar(value.(phi2_sh[ens_340]), value.(t0fk_st_sh[ens_340]), err.(t0fk_st_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(t0fk_st_sh[ens_346]), err.(t0fk_st_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(t0fk_st_sh[ens_355]), err.(t0fk_st_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(t0fk_st_sh[ens_370]), err.(t0fk_st_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(t0fk_st_sh[ens_385]), err.(t0fk_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    ax[:set_ylim]([0.3, 0.34])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0fk_st.pdf")

    #t0fpi_sh = sqrt.(t0_sh) .* fpi_extrap
    uwerr.(t0fpi_sh)
    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{8t_0}f_{\pi}$")
    errorbar(value.(phi2_sh[ens_340]), value.(t0fpi_sh[ens_340]), err.(t0fpi_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(t0fpi_sh[ens_346]), err.(t0fpi_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(t0fpi_sh[ens_355]), err.(t0fpi_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(t0fpi_sh[ens_370]), err.(t0fpi_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(t0fpi_sh[ens_385]), err.(t0fpi_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    ax[:set_ylim]([0.24, 0.321])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0fpi_tm.pdf")

    #t0fk_sh = sqrt.(t0_sh) .* fk_extrap
    uwerr.(t0fk_sh)
    uwerr.(phi2_sh)
    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{8t_0}f_{K}$")
    errorbar(value.(phi2_sh[ens_340]), value.(t0fk_sh[ens_340]), err.(t0fk_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346]), value.(t0fk_sh[ens_346]), err.(t0fk_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355]), value.(t0fk_sh[ens_355]), err.(t0fk_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370]), value.(t0fk_sh[ens_370]), err.(t0fk_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385]), value.(t0fk_sh[ens_385]), err.(t0fk_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    ax = gca()
    ax[:set_ylim]([0.3, 0.34])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0fk_tm.pdf")

#==============================================================================================================================#

#============================== continuum limit symmetric point ===============================================================#

    function model_continuum(x,p) 
        return [p[1] + p[2] * x[i] for i in 1:length(x)]
    end

    x = 1 ./ (8 * t0_sh[ind_sym])
    uwerr.(x)
    y = t0fpik_sh[ind_sym]
    y_st = t0fpik_st_sh[ind_sym]

    up, chi2, chi_exp, pv = fit_alg(model_continuum,value.(x),y,2)
    up_st, chi2_st, chi_exp_st, pv_st = fit_alg(model_continuum,value.(x),y_st,2)

    aux = model_continuum(x, up) ; uwerr.(aux)
    aux_st = model_continuum(x, up_st) ; uwerr.(aux_st)

    fig = figure("pyplot_subplot_column10")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    errorbar(value(x[1]), value(y[1]), err(y[1]), err(x[1]), fmt="s", label=L"$\beta=3.40$", color="rebeccapurple")
    errorbar(value(x[2]), value(y[2]), err(y[2]), err(x[2]), fmt="o",  label=L"$\beta=3.46$", color="green")
    errorbar(value(x[3]), value(y[3]), err(y[3]), err(x[3]), fmt="<",  label=L"$\beta=3.55$", color="blue")
    errorbar(value(x[4]), value(y[4]), err(y[4]), err(x[4]), fmt=">",  label=L"$\beta=3.70$", color="darkorange")
    errorbar(value(x[5]), value(y[5]), err(y[5]), err(x[5]), fmt=">",  label=L"$\beta=3.85$", color="red")
    errorbar(value(x[1]), value(y_st[1]), err(y_st[1]), err(x[1]), fmt="s", mfc="none", color="rebeccapurple")
    errorbar(value(x[2]), value(y_st[2]), err(y_st[2]), err(x[2]), fmt="o", mfc="none", color="green")
    errorbar(value(x[3]), value(y_st[3]), err(y_st[3]), err(x[3]), fmt="<", mfc="none", color="blue")
    errorbar(value(x[4]), value(y_st[4]), err(y_st[4]), err(x[4]), fmt=">", mfc="none", color="darkorange")
    errorbar(value(x[5]), value(y_st[5]), err(y_st[5]), err(x[5]), fmt=">", mfc="none", color="red")
    xlabel(L"$a^2/8t_0$")
    ylabel(L"$\sqrt{8t_0}f_{\pi K}$")
    x_plot = [i for i in 0.00:0.0005:0.045]
    aux = model_continuum(x_plot, up) ; uwerr.(aux)
    aux_st = model_continuum(x_plot, up_st) ; uwerr.(aux_st)
    v = value.(aux)
    e = err.(aux)
    fill_between(x_plot, v-e, v+e, color="fuchsia", alpha=0.1)
    v_2 = value.(aux_st)
    e_2 = err.(aux_st)
    fill_between(x_plot, v_2-e_2, v_2+e_2, color="fuchsia", alpha=0.2)
    ax = gca()
    #ax[:set_ylim]([0.295, 0.325])
    legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/continuum_sym.pdf")

#==============================================================================================================================#

#============================== chiral & continuum limits =====================================================================#

    ## prepare data SU3 & Tay
        TIC = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        chi = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        W = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        t0fpik_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        pval = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        uprm_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_plot2 = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_plot_SU2 = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]

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

        y = t0fpik_sh
        y_st = t0fpik_st_sh
        y_combined = [t0fpik_st_sh; t0fpik_sh]
        y_355 = y[[ens_346; ens_355; ens_370; ens_385]]
        y_st_355 = y_st[[ens_346; ens_355; ens_370; ens_385]]
        y_combined_355 = [y_st_355; y_355]
        y_3552 = y[[ens_355; ens_370; ens_385]]
        y_st_3552 = y_st[[ens_355; ens_370; ens_385]]
        y_combined_3552 = [y_st_3552; y_3552]
        y_combined_355_only_Wil = [y_st_355; y]
        y_combined_355_only_Wtm = [y_st; y_355]
        y_combined_355 = [y_st_355; y_355]
        y_nosym = y[ind_nosym]
        y_st_nosym = y_st[ind_nosym]
        y_combined_nosym = [y_st_nosym; y_nosym]
        y_nosym355 = y[ind_nosym355]
        y_st_nosym355 = y_st[ind_nosym355]
        y_combined_nosym355 = [y_st_nosym355; y_nosym355]
        y_phi204 = y[ind_phi204]
        y_st_phi204 = y_st[ind_phi204]
        y_combined_phi204 = [y_st_phi204; y_phi204]
        y_mL_40 = y[ind_mL_40]
        y_mL_41 = y[ind_mL_41]
        y_mL_42 = y[ind_mL_42]
        y_st_mL_40 = y_st[ind_mL_40]
        y_st_mL_41 = y_st[ind_mL_41]
        y_st_mL_42 = y_st[ind_mL_42]
        y_combined_mL_40 = [y_st_mL_40; y_mL_40]
        y_combined_mL_41 = [y_st_mL_41; y_mL_41]
        y_combined_mL_42 = [y_st_mL_42; y_mL_42]

        cuts_y = [y, y_355, y_3552, y_nosym, y_nosym355, y_phi204, y_mL_42]
        cuts_y_st = [y_st, y_st_355, y_st_3552, y_st_nosym, y_st_nosym355, y_st_phi204, y_st_mL_42]
        cuts_y_combined = [y_combined, y_combined_355, y_combined_3552, y_combined_nosym, y_combined_nosym355, y_combined_phi204, y_combined_mL_42]
        cuts_x = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_st = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_combined = [x_combined, x_combined_355, x_combined_3552, x_combined_nosym, x_combined_nosym355, x_combined_phi204, x_combined_mL_42]

        set_y = [cuts_y, cuts_y_st, cuts_y_combined]
        set_x = [cuts_x, cuts_x_st, cuts_x_combined]

        #Wm = [inv.(Symmetric.(cov.(set_y[i]))) for i in 1:length(set_y)]
    ##

    switch = 0 ## .02-2
    #Wm_syst = [inv.(Symmetric.(diagm.(diag.(cov.(set_y[i]))) .+ switch ^ 2 * [diagm(value.(set_x[i][k][:,1] .^ 4)) for k in 1:length(set_y[i])])) for i in 1:length(set_y)]    
    #Wm_syst = [convert.(Matrix{Float64}, Wm_syst[i]) for i in 1:length(set_y)]
    #Wm = [inv.(Symmetric.(diagm.(diag.(cov.(set_y[i]))))) for i in 1:length(set_y)]
    #Wm = [convert.(Matrix{Float64}, Wm[i]) for i in 1:length(set_y)]
    Wm_syst = [inv.(Symmetric.(((cov.(set_y[i]))) .+ switch ^ 2 * [diagm(value.(set_x[i][k][:,1] .^ 4)) for k in 1:length(set_y[i])])) for i in 1:length(set_y)]    
    Wm_syst = [convert.(Matrix{Float64}, Wm_syst[i]) for i in 1:length(set_y)]
    Wm = [inv.(Symmetric.(((cov.(set_y[i]))))) for i in 1:length(set_y)]
    Wm = [convert.(Matrix{Float64}, Wm[i]) for i in 1:length(set_y)]

    models = [model2_ChPT_a2; model2_ChPT_aas; model2_ChPT_a2phi2; model2_Taylor_a2; model2_Taylor_aas; model2_Taylor_a2phi2; model2_Taylor4_a2; model2_Taylor4_a2phi2]
    models_combined = [model2_ChPT_a2_combined; model2_ChPT_aas_combined; model2_ChPT_a2a2phi2_combined; model2_ChPT_a2phi2a2_combined; model2_ChPT_a2phi2_combined; model2_Taylor_a2_combined; model2_Taylor_aas_combined; model2_Taylor_a2a2phi2_combined; model2_Taylor_a2phi2a2_combined; model2_Taylor_a2phi2_combined; model2_Taylor4_a2_combined; model2_Taylor4_a2a2phi2_combined; model2_Taylor4_a2phi2a2_combined; model2_Taylor4_a2phi2_combined]
    models = [models, models, models_combined]
    param = [3,3,4,3,3,4,4,5]
    param_combined = [4,4,5,5,6,4,4,5,5,6,5,6,6,7]
    param = [param, param, param_combined]

    for k in 1:length(set_y)
        for i in 1:length(models[k])
            for j in 1:length(cuts_y)
                if (k in [1,2] && j in [5,6] && i in [3,6,8] || k == 3 && j in [5,6] && i in [3,4,5,8,9,10,12,13,14]) == false
                    x = set_x[k][j]
                    y = set_y[k][j]
                    global L1 = length(set_y[1][j])
                    global L2 = length(set_y[1][j])
                    uprm, chi2, chi_exp, pval_aux, doff = fit_alg(models[k][i], value.(x), y, param[k][i], Wm_syst[k][j])
                    push!(TIC[k], chi2 - 2 * chi_exp)
                    push!(pval[k], pval_aux)
                    push!(t0fpik_ph_vec[k], models[k][i]([x_ph;x],uprm)[1])
                    if i == j == 1 uprm_plot[k] = uprm end
                    if i == 3 && j == 1 uprm_plot2[k] = uprm end
                end
            end
        end
    end

    ## prepare data SU2 only fpi
        uprm_plot_SU2 = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        t0fpi_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]

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
        x_mL_41 = x[ind_mL_41, :]
        x_mL_42 = x[ind_mL_42, :]
        x_combined_mL_41 = [x_mL_41; x_mL_41]
        x_combined_mL_42 = [x_mL_42; x_mL_42]
        x_ph = [0.0 (phi2_ph) (phi4_ph) (phi2_sym_ph)]
        x_3552 = x[[ens_355; ens_370; ens_385],:]
        x_combined_3552 = [x_3552; x_3552]
        x_phi204 = x[ind_phi204,:]
        x_combined_phi204 = [x_phi204; x_phi204]

        fpi = t0fpi_sh
        fpi_355 = t0fpi_sh[[ens_346; ens_355; ens_370; ens_385]]
        fpi_3552 = t0fpi_sh[[ens_355; ens_370; ens_385]]
        fpi_nosym = t0fpi_sh[ind_nosym]
        fpi_nosym355 = t0fpi_sh[ind_nosym355]
        fpi_phi204 = t0fpi_sh[ind_phi204]
        fpi_mL_41 = t0fpi_sh[ind_mL_41]
        fpi_mL_42 = t0fpi_sh[ind_mL_42]
        fpi_st = t0fpi_st_sh
        fpi_st_355 = t0fpi_st_sh[[ens_346; ens_355; ens_370; ens_385]]
        fpi_st_3552 = t0fpi_st_sh[[ens_355; ens_370; ens_385]]
        fpi_st_nosym = t0fpi_st_sh[ind_nosym]
        fpi_st_nosym355 = t0fpi_st_sh[ind_nosym355]
        fpi_st_phi204 = t0fpi_st_sh[ind_phi204]
        fpi_st_mL_41 = t0fpi_st_sh[ind_mL_41]
        fpi_st_mL_42 = t0fpi_st_sh[ind_mL_42]
        y = fpi
        y_st = fpi_st
        y_combined = [y_st; y]
        y_355 = fpi_355
        y_st_355 = fpi_st_355
        y_combined_355 = [y_st_355; y_355]
        y_3552 = fpi_3552
        y_st_3552 = fpi_st_3552
        y_combined_3552 = [y_st_3552; y_3552]
        y_nosym = fpi_nosym
        y_nosym355 = fpi_nosym355
        y_st_nosym = fpi_st_nosym
        y_st_nosym355 = fpi_st_nosym355
        y_combined_nosym = [y_st_nosym; y_nosym]
        y_combined_nosym355 = [y_st_nosym355; y_nosym355]
        y_phi204 = fpi_phi204
        y_st_phi204 = fpi_st_phi204
        y_combined_phi204 = [y_st_phi204; y_phi204]
        y_mL_41 = fpi_mL_41
        y_st_mL_41 = fpi_st_mL_41
        y_combined_mL_41 = [y_st_mL_41; y_mL_41]
        y_mL_42 = fpi_mL_42
        y_st_mL_42 = fpi_st_mL_42
        y_combined_mL_42 = [y_st_mL_42; y_mL_42]

        cuts_y = [y, y_355, y_3552, y_nosym, y_nosym355, y_phi204, y_mL_42]
        cuts_y_st = [y_st, y_st_355, y_st_3552, y_st_nosym, y_st_nosym355, y_st_phi204, y_st_mL_42]
        cuts_y_combined = [y_combined, y_combined_355, y_combined_3552, y_combined_nosym, y_combined_nosym355, y_combined_phi204, y_combined_mL_42]
        
        cuts_x = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_st = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_combined = [x_combined, x_combined_355, x_combined_3552, x_combined_nosym, x_combined_nosym355, x_combined_phi204, x_combined_mL_42]

        set_y = [cuts_y, cuts_y_st, cuts_y_combined]
        set_x = [cuts_x, cuts_x_st, cuts_x_combined]
        [[uwerr.(set_y[i][j]) for j in 1:length(set_y[i])] for i in 1:length(set_y)]
        [[uwerr.(set_x[i][j]) for j in 1:length(set_x[i])] for i in 1:length(set_y)]
    ##

    switch = 0 ## .02-2
    #Wm_syst = [inv.(Symmetric.(diagm.(diag.(cov.(set_y[i]))) .+ switch ^ 2 * [diagm(value.(set_x[i][k][:,1] .^ 4)) for k in 1:length(set_y[i])])) for i in 1:length(set_y)]    
    #Wm_syst = [convert.(Matrix{Float64}, Wm_syst[i]) for i in 1:length(set_y)]
    #Wm = [inv.(Symmetric.(diagm.(diag.(cov.(set_y[i]))))) for i in 1:length(set_y)]
    #Wm = [convert.(Matrix{Float64}, Wm[i]) for i in 1:length(set_y)]
    Wm_syst = [inv.(Symmetric.(((cov.(set_y[i]))) .+ switch ^ 2 * [diagm(value.(set_x[i][k][:,1] .^ 4)) for k in 1:length(set_y[i])])) for i in 1:length(set_y)]    
    Wm_syst = [convert.(Matrix{Float64}, Wm_syst[i]) for i in 1:length(set_y)]
    Wm = [inv.(Symmetric.(((cov.(set_y[i]))))) for i in 1:length(set_y)]
    Wm = [convert.(Matrix{Float64}, Wm[i]) for i in 1:length(set_y)]

    models = [model_ChPT2_a2; model_ChPT2_aas; model_ChPT2_a2phi2]
    models_combined = [model_ChPT2_a2_combined; model_ChPT2_aas_combined; model_ChPT2_a2a2phi2_combined; model_ChPT2_a2phi2a2_combined; model_ChPT2_a2_combined]
    models = [models, models, models_combined]
    param = [4,4,5]
    param_combined = [5,5,6,6,7]
    param = [param, param, param_combined]

    for k in 1:length(set_y)
        for i in 1:length(models[k])
            for j in 1:length(cuts_y)
                if (k in [1,2] && j in [5,6] && i in [2] || k == 3 && j in [5,6] && i in [2,3,4]) == false
                    x = set_x[k][j]
                    y = set_y[k][j]
                    global L1 = length(set_y[1][j])
                    global L2 = length(set_y[1][j])
                    uprm, chi2, chi_exp, pval_aux, doff = fit_alg(models[k][i], value.(x), y, param[k][i], Wm_syst[k][j])
                    push!(TIC[k], chi2 - 2 * chi_exp)
                    push!(pval[k], pval_aux)
                    if k == 3
                        push!(t0fpi_ph_vec[k], model_plot_SU2_pi([x_ph;x_ph;x_ph;x_ph],uprm)[1])
                    else
                        push!(t0fpi_ph_vec[k], model_plot_SU2_pi([x_ph;x_ph],uprm)[1])
                    end
                    if i == 1 && j == 1 uprm_plot_SU2[k] = uprm end
                end
            end
        end
    end

    ## prepare data SU2 fpi & fk
        uprm_plot_SU2 = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        t0fk_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        t0fpi_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]

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
        x_mL_41 = x[ind_mL_41, :]
        x_mL_42 = x[ind_mL_42, :]
        x_combined_mL_41 = [x_mL_41; x_mL_41]
        x_combined_mL_42 = [x_mL_42; x_mL_42]
        x_ph = [0.0 (phi2_ph) (phi4_ph) (phi2_sym_ph)]
        x_3552 = x[[ens_355; ens_370; ens_385],:]
        x_combined_3552 = [x_3552; x_3552]
        x_phi204 = x[ind_phi204,:]
        x_combined_phi204 = [x_phi204; x_phi204]

        fpi = t0fpi_sh
        fpi_355 = t0fpi_sh[[ens_346; ens_355; ens_370; ens_385]]
        fpi_3552 = t0fpi_sh[[ens_355; ens_370; ens_385]]
        fpi_nosym = t0fpi_sh[ind_nosym]
        fpi_nosym355 = t0fpi_sh[ind_nosym355]
        fpi_phi204 = t0fpi_sh[ind_phi204]
        fpi_mL_41 = t0fpi_sh[ind_mL_41]
        fpi_mL_42 = t0fpi_sh[ind_mL_42]
        fk = t0fk_sh
        fk_355 = t0fk_sh[[ens_346; ens_355; ens_370; ens_385]]
        fk_3552 = t0fk_sh[[ens_355; ens_370; ens_385]]
        fk_nosym = t0fk_sh[ind_nosym]
        fk_nosym355 = t0fk_sh[ind_nosym355]
        fk_phi204 = t0fk_sh[ind_phi204]
        fk_mL_41 = t0fk_sh[ind_mL_41]
        fk_mL_42 = t0fk_sh[ind_mL_42]
        fpi_st = t0fpi_st_sh
        fpi_st_355 = t0fpi_st_sh[[ens_346; ens_355; ens_370; ens_385]]
        fpi_st_3552 = t0fpi_st_sh[[ens_355; ens_370; ens_385]]
        fpi_st_nosym = t0fpi_st_sh[ind_nosym]
        fpi_st_nosym355 = t0fpi_st_sh[ind_nosym355]
        fpi_st_phi204 = t0fpi_st_sh[ind_phi204]
        fpi_st_mL_41 = t0fpi_st_sh[ind_mL_41]
        fpi_st_mL_42 = t0fpi_st_sh[ind_mL_42]
        fk_st = t0fk_st_sh
        fk_st_355 = t0fk_st_sh[[ens_346; ens_355; ens_370; ens_385]]
        fk_st_3552 = t0fk_st_sh[[ens_355; ens_370; ens_385]]
        fk_st_nosym = t0fk_st_sh[ind_nosym]
        fk_st_nosym355 = t0fk_st_sh[ind_nosym355]
        fk_st_phi204 = t0fk_st_sh[ind_phi204]
        fk_st_mL_41 = t0fk_st_sh[ind_mL_41]
        fk_st_mL_42 = t0fk_st_sh[ind_mL_42]
        y = [fpi; fk]
        y_st = [fpi_st; fk_st]
        y_combined = [y_st; y]
        y_355 = [fpi_355; fk_355]
        y_st_355 = [fpi_st_355; fk_st_355]
        y_combined_355 = [y_st_355; y_355]
        y_3552 = [fpi_3552; fk_3552]
        y_st_3552 = [fpi_st_3552; fk_st_3552]
        y_combined_3552 = [y_st_3552; y_3552]
        y_nosym = [fpi_nosym; fk_nosym]
        y_nosym355 = [fpi_nosym355; fk_nosym355]
        y_st_nosym = [fpi_st_nosym; fk_st_nosym]
        y_st_nosym355 = [fpi_st_nosym355; fk_st_nosym355]
        y_combined_nosym = [y_st_nosym; y_nosym]
        y_combined_nosym355 = [y_st_nosym355; y_nosym355]
        y_phi204 = [fpi_phi204; fk_phi204]
        y_st_phi204 = [fpi_st_phi204; fk_st_phi204]
        y_combined_phi204 = [y_st_phi204; y_phi204]
        y_mL_41 = [fpi_mL_41; fk_mL_41]
        y_st_mL_41 = [fpi_st_mL_41; fk_st_mL_41]
        y_combined_mL_41 = [y_st_mL_41; y_mL_41]
        y_mL_42 = [fpi_mL_42; fk_mL_42]
        y_st_mL_42 = [fpi_st_mL_42; fk_st_mL_42]
        y_combined_mL_42 = [y_st_mL_42; y_mL_42]

        cuts_y = [y, y_355, y_3552, y_nosym, y_nosym355, y_phi204, y_mL_42]

        cuts_y = [y, y_355, y_3552, y_nosym, y_nosym355, y_phi204, y_mL_42]
        cuts_y_st = [y_st, y_st_355, y_st_3552, y_st_nosym, y_st_nosym355, y_st_phi204, y_st_mL_42]
        cuts_y_combined = [y_combined, y_combined_355, y_combined_3552, y_combined_nosym, y_combined_nosym355, y_combined_phi204, y_combined_mL_42]
        
        cuts_x = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_st = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_combined = [x_combined, x_combined_355, x_combined_3552, x_combined_nosym, x_combined_nosym355, x_combined_phi204, x_combined_mL_42]
        cuts_x = [[cuts_x[i]; cuts_x[i]] for i in 1:length(cuts_x)]
        cuts_x_st = [[cuts_x_st[i]; cuts_x_st[i]] for i in 1:length(cuts_x_st)]
        cuts_x_combined = [[cuts_x_combined[i]; cuts_x_combined[i]] for i in 1:length(cuts_x_combined)]

        set_y = [cuts_y, cuts_y_st, cuts_y_combined]
        set_x = [cuts_x, cuts_x_st, cuts_x_combined]
        [[uwerr.(set_y[i][j]) for j in 1:length(set_y[i])] for i in 1:length(set_y)]
        [[uwerr.(set_x[i][j]) for j in 1:length(set_x[i])] for i in 1:length(set_y)]
    ##

    models = [model2_ChPT2_a2; model2_ChPT2_aas; model2_ChPT2_a2phi2]
    models_combined = [model2_ChPT2_a2_combined; model2_ChPT2_aas_combined; model2_ChPT2_a2a2phi2_combined; model2_ChPT2_a2phi2a2_combined; model2_ChPT2_a2_combined]
    models = [models, models, models_combined]
    param = [7,7,9]
    param_combined = [9,9,11,11,13]
    param = [param, param, param_combined]

    for k in 1:length(set_y)
        for i in 1:length(models[k])
            for j in 1:length(cuts_y)
                if (k in [1,2] && j in [5,6] && i in [2] || k == 3 && j in [5,6] && i in [2,3,4]) == false
                    x = set_x[k][j]
                    y = set_y[k][j]
                    global L1 = length(set_y[1][j])
                    global L2 = length(set_y[1][j])
                    guess = [-0.03714930996391645, 3.794500859421803, -0.04370336668443585, -0.3767099865988463, -0.053522736192504715, 0.32081397340658246, -0.340720745616767]
                    uprm, chi2, chi_exp, pval_aux = fit_alg(models[k][i], value.(x), y, param[k][i], guess)
                    push!(TIC[k], chi2 - 2 * chi_exp)
                    push!(pval[k], pval_aux)
                    if k == 3
                        push!(t0fpi_ph_vec[k], model_plot_SU2_pi([x_ph;x_ph;x_ph;x_ph],uprm)[1])
                        push!(t0fk_ph_vec[k], model_plot_SU2_k([x_ph;x_ph;x_ph;x_ph],uprm)[2])
                    else
                        push!(t0fpi_ph_vec[k], model_plot_SU2_pi([x_ph;x_ph],uprm)[1])
                        push!(t0fk_ph_vec[k], model_plot_SU2_k([x_ph;x_ph],uprm)[2])
                    end
                    if i == 1 && j == 1 uprm_plot_SU2[k] = uprm end
                end
            end
        end
    end

    TIC = [TIC[k] .- minimum.(TIC)[k] for k in 1:length(TIC)]
    W = [exp.(-0.5 * TIC[k]) ./ sum(exp.(-0.5 * TIC[k])) for k in 1:length(TIC)]
    #t0fpik_ph_vec = [[t0fpik_ph_vec[k]; 2/3 * (t0fk_ph_vec .+ 0.5 * t0fpi_ph_vec)[k]] for k in 1:length(t0fpik_ph_vec)]
    t0fpik_ph = [sum(t0fpik_ph_vec[k] .* W[k]) for k in 1:length(TIC)]
    syst = [sqrt(sum(t0fpik_ph_vec[k] .^ 2 .* W[k]) - (sum(t0fpik_ph_vec[k] .* W[k])) ^ 2) for k in 1:length(TIC)]
    t0fpik_ph = t0fpik_ph .+ [uwreal([0.0, value(syst[k])], string("syst chiral ", k, " 3rd")) for k in 1:length(TIC)]
    uwerr.(t0fpik_ph)
    sqrt_t0_ph = [t0fpik_ph[k] / (sqrt(8) * fpik_exp) for k in 1:length(TIC)]
    uwerr.(sqrt_t0_ph)

    #=
    W_fpi = [W[k][end-length(t0fpi_ph_vec[k])+1:end] for k in 1:3]
    pval_fpi = [pval[k][end-length(t0fpi_ph_vec[k])+1:end] for k in 1:3]
    W_fpi = [W_fpi[k] ./ sum(W_fpi[k]) for k in 1:3]
    t0fpi_ph = [sum(t0fpi_ph_vec[k] .* W_fpi[k]) for k in 1:length(TIC)]
    syst_fpi = [sqrt(sum(t0fpi_ph_vec[k] .^ 2 .* W_fpi[k]) - (sum(t0fpi_ph_vec[k] .* W_fpi[k])) ^ 2) for k in 1:length(TIC)]
    t0fpi_ph = t0fpi_ph .+ [uwreal([0.0, value(syst_fpi[k])], string("syst chiral fpi", k, " 3rd")) for k in 1:length(TIC)]
    uwerr.(t0fpi_ph)
    sqrt_t0_ph_fpi = [t0fpi_ph[k] / (sqrt(8) * Fpi / hc) for k in 1:length(TIC)]
    uwerr.(sqrt_t0_ph_fpi)
    t0fk_ph = [sum(t0fk_ph_vec[k] .* W_fpi[k]) for k in 1:length(TIC)]
    syst_fk = [sqrt(sum(t0fk_ph_vec[k] .^ 2 .* W_fpi[k]) - (sum(t0fk_ph_vec[k] .* W_fpi[k])) ^ 2) for k in 1:length(TIC)]
    t0fk_ph = t0fk_ph .+ [uwreal([0.0, value(syst_fk[k])], string("syst chiral fpi ", k, " 3rd")) for k in 1:length(TIC)]
    fk_exp_pred = t0fk_ph ./ sqrt(8) ./ sqrt_t0_ph_fpi * hc
    uwerr.(fk_exp_pred)
    W_fpi_aux = deepcopy(W_fpi)
    pval_fpi_aux = deepcopy(pval_fpi)
    [uwerr.(t0fpi_ph_vec[k]) for k in 1:length(TIC)]
    [uwerr.(t0fk_ph_vec[k]) for k in 1:length(TIC)]
    sqrt_t0_ph_fpi_vec = [[t0fpi_ph_vec[k][i] / (sqrt(8)* Fpi / hc) for i in 1:length(t0fpi_ph_vec[k])] for k in 1:length(W_fpi)]
    [uwerr.(sqrt_t0_ph_fpi_vec[k]) for k in 1:length(TIC)]
    =#

    W_aux = deepcopy(W)
    pval_aux = deepcopy(pval)
    [uwerr.(t0fpik_ph_vec[k]) for k in 1:length(TIC)]
    [uwerr.(t0fpi_ph_vec[k]) for k in 1:length(TIC)]
    [uwerr.(t0fk_ph_vec[k]) for k in 1:length(TIC)]

    sqrt_t0_ph_vec = [[t0fpik_ph_vec[k][i] / (sqrt(8)* fpik_exp) for i in 1:length(t0fpik_ph_vec[k])] for k in 1:length(TIC)]
    [uwerr.(sqrt_t0_ph_vec[k]) for k in 1:length(TIC)]
    ixx = 3
    details(sqrt_t0_ph[ixx], string("syst chiral ",ixx," 3rd")); sqrt(1-54/100) * err(sqrt_t0_ph[ixx])

    ## à la Strassberger:
    ixx = 3
    ix = findall(pval[ixx] .> .1)
    a = (maximum(value.(sqrt_t0_ph_vec[ixx][ix])) - minimum(value.(sqrt_t0_ph_vec[ixx][ix]))) / 2
    b = maximum(pval[ixx])
    c = findall(pval[ixx] .>= b)
    sqrt_t0_ph_vec[ixx][c]

    #=
    fb = BDIO_open("/home/asaez/cls_ens/results/t0_combined_4th.bdio", "w")
    write_uwreal(sqrt_t0_ph[1] ^ 2, fb, 1)
    write_uwreal(sqrt_t0_ph[2] ^ 2, fb, 2)
    write_uwreal(sqrt_t0_ph[3] ^ 2, fb, 3)
    write_uwreal(sqrt(8) * sqrt_t0_ph[3] * fpik_exp, fb, 4)
    BDIO_close!(fb)
    =#

    phi4_new = 8 * sqrt_t0_ph[3] ^ 2 * (mk_exp ^ 2 + 0.5 * mpi_exp ^ 2)

#==============================================================================================================================#

#============================== plot chiral & continuum =======================================================================#

    uwerr(phi2_ph)
    [uwerr.(t0fpik_ph_vec[k]) for k in 1:length(TIC)]
    [uwerr.(t0fpi_ph_vec[k]) for k in 1:length(TIC)]
    [uwerr.(t0fk_ph_vec[k]) for k in 1:length(TIC)]

    #t0fpik SU3 only combined
        uprm_combined = uprm_plot[3]
        uprm_st = uprm_plot[2]
        uprm = uprm_plot[1]

        Fph = model_plot([0 ./ t0_sh ./ 8 [phi2_ph for i in 1:length(t0_sh)] phi4_sh], uprm_combined[[1,2,4]])
        Fphi2 = model_plot([0 ./ t0_sh ./ 8 phi2_sh phi4_sh], uprm_combined[[1,2,4]])
        y_aux = t0fpik_sh .- (Fphi2 .- Fph)
        Fph_st = model_plot([0 ./ t0_sh ./ 8 [phi2_ph for i in 1:length(t0_sh)] phi4_sh], uprm_combined[[1,2,3]])
        Fphi2_st = model_plot([0 ./ t0_sh ./ 8 phi2_sh phi4_sh], uprm_combined[[1,2,3]])
        y_aux_st = t0fpik_st_sh .- (Fphi2_st .- Fph_st)
        uwerr.(y_aux)
        a2t0 = 1 ./ t0_sh ./ 8
        uwerr.(a2t0)

        fig = figure("SU3")
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 15

        color_beta = ["rebeccapurple", "green", "blue", "darkorange", "red"]

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{8t_0}f_{\pi K}$")
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_st_sh[ens_340]), err.(t0fpik_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_st_sh[ens_346]), err.(t0fpik_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_st_sh[ens_355]), err.(t0fpik_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_st_sh[ens_370]), err.(t0fpik_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_st_sh[ens_385]), err.(t0fpik_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_sh[ens_340]), err.(t0fpik_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_sh[ens_346]), err.(t0fpik_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_sh[ens_355]), err.(t0fpik_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_sh[ens_370]), err.(t0fpik_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_sh[ens_385]), err.(t0fpik_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            list = [1,2,4]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            #fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot(x_plot,uprm_combined) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            list = [1,2,3]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpik_ph_vec[3][1]), err(t0fpik_ph_vec[3][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.283, 0.325])

        subplot(122)
        errorbar(value.(a2t0[ens_340]), value.(y_aux[ens_340]), err.(y_aux[ens_340]), err.(a2t0[ens_340]), fmt="s", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(a2t0[ens_346]), value.(y_aux[ens_346]), err.(y_aux[ens_346]), err.(a2t0[ens_346]), fmt="o", label=L"$\beta=3.46$", color="green")
        errorbar(value.(a2t0[ens_355]), value.(y_aux[ens_355]), err.(y_aux[ens_355]), err.(a2t0[ens_355]), fmt="<", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(a2t0[ens_370]), value.(y_aux[ens_370]), err.(y_aux[ens_370]), err.(a2t0[ens_370]), fmt=">", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(a2t0[ens_385]), value.(y_aux[ens_385]), err.(y_aux[ens_385]), err.(a2t0[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$a^2/8t_0$")
        x_prime = [i for i in 0.0:0.001:0.05]
        x_plot = [x_prime [phi2_ph for i in 1:length(x_prime)] [(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot(x_plot,uprm_combined[[1,2,4]]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.75)
        uwerr.(y_aux_st)
        a2t0 = 1 ./ t0_sh ./ 8
        uwerr.(a2t0)
        errorbar(value.(a2t0[ens_340]), value.(y_aux_st[ens_340]), err.(y_aux_st[ens_340]), err.(a2t0[ens_340]), fmt="s", mfc="none", color="rebeccapurple")
        errorbar(value.(a2t0[ens_346]), value.(y_aux_st[ens_346]), err.(y_aux_st[ens_346]), err.(a2t0[ens_346]), fmt="o", mfc="none", color="green")
        errorbar(value.(a2t0[ens_355]), value.(y_aux_st[ens_355]), err.(y_aux_st[ens_355]), err.(a2t0[ens_355]), fmt="<", mfc="none", color="blue")
        errorbar(value.(a2t0[ens_370]), value.(y_aux_st[ens_370]), err.(y_aux_st[ens_370]), err.(a2t0[ens_370]), fmt=">", mfc="none", color="darkorange")
        errorbar(value.(a2t0[ens_385]), value.(y_aux_st[ens_385]), err.(y_aux_st[ens_385]), err.(a2t0[ens_385]), fmt="^", mfc="none", color="red")
        xlabel(L"$a^2/8t_0$")
        x_prime = [i for i in 0.0:0.001:0.05]
        x_plot = [x_prime [phi2_ph for i in 1:length(x_prime)] [(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot(x_plot,uprm_combined[[1,2,3]]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)

        errorbar(0, value(t0fpik_ph_vec[3][1]), err(t0fpik_ph_vec[3][1]), 0, fmt="x", label="ph. point", color="black")
        errorbar(-0.001, value(t0fpik_ph[3]), err(t0fpik_ph[3]), 0, fmt="*", label="ph. point model av", color="black")
        ax = gca()
        ax[:set_ylim]([0.283, 0.325])
        setp(ax.get_yticklabels(),visible=false)
        #legend(loc="upper right", bbox_to_anchor=(3.,1.))
        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/SU3_comb.pdf")

    ##

    #t0fpik SU3
        uprm_combined = uprm_plot[3]
        uprm_st = uprm_plot[2]
        uprm = uprm_plot[1]

        fig = figure("SU3")
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        color_beta = ["rebeccapurple", "green", "blue", "darkorange", "red"]

        subplot(131) # Create the 1st axis of a 3x1 array of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_sh[ens_340]), err.(t0fpik_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_sh[ens_346]), err.(t0fpik_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_sh[ens_355]), err.(t0fpik_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_sh[ens_370]), err.(t0fpik_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_sh[ens_385]), err.(t0fpik_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{8t_0}f_{\pi K}$")
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot(x_plot,uprm) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot(x_plot,uprm) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpik_ph_vec[1][1]), err(t0fpik_ph_vec[1][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.283, 0.325])
        legend()
        tight_layout()
            
        subplot(132) # Create the 2nd axis of a 3x1 arrax of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_st_sh[ens_340]), err.(t0fpik_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_st_sh[ens_346]), err.(t0fpik_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_st_sh[ens_355]), err.(t0fpik_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_st_sh[ens_370]), err.(t0fpik_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_st_sh[ens_385]), err.(t0fpik_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot(x_plot,uprm_st) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot(x_plot,uprm_st) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpik_ph_vec[2][1]), err(t0fpik_ph_vec[2][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.283, 0.325])
        setp(ax.get_yticklabels(),visible=false)

        subplot(133) # Create the 3rd axis of a 3x1 array of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_st_sh[ens_340]), err.(t0fpik_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_st_sh[ens_346]), err.(t0fpik_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_st_sh[ens_355]), err.(t0fpik_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_st_sh[ens_370]), err.(t0fpik_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_st_sh[ens_385]), err.(t0fpik_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpik_sh[ens_340]), err.(t0fpik_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpik_sh[ens_346]), err.(t0fpik_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpik_sh[ens_355]), err.(t0fpik_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpik_sh[ens_370]), err.(t0fpik_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpik_sh[ens_385]), err.(t0fpik_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        i = 1
        for ind in ind_sym
            list = [1,2,4]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot(x_plot,uprm_combined) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            list = [1,2,3]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpik_ph_vec[3][1]), err(t0fpik_ph_vec[3][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.283, 0.325])
        setp(ax.get_yticklabels(),visible=false)
        subplots_adjust(wspace=0.02)

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/SU3.pdf")
    ##

    #t0fpi SU2
        uprm_combined = uprm_plot_SU2[3]
        uprm_st = uprm_plot_SU2[2]
        uprm = uprm_plot_SU2[1]

        fig = figure("SU2 pion")
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        color_beta = ["rebeccapurple", "green", "blue", "darkorange", "red"]

        subplot(131) # Create the 1st axis of a 3x1 array of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpi_sh[ens_340]), err.(t0fpi_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpi_sh[ens_346]), err.(t0fpi_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpi_sh[ens_355]), err.(t0fpi_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpi_sh[ens_370]), err.(t0fpi_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpi_sh[ens_385]), err.(t0fpi_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{8t_0}f_{\pi}$")
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot_SU2_pi(x_plot,uprm) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_pi(x_plot,uprm) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            #fill_between(x_plot[:,2], v-e, v+e, color_beta[i], alpha=0.75)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpi_ph_vec[1][1]), err(t0fpi_ph_vec[1][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.24, 0.33])
        #legend()
        tight_layout()
            
        subplot(132) # Create the 2nd axis of a 3x1 arrax of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpi_st_sh[ens_340]), err.(t0fpi_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpi_st_sh[ens_346]), err.(t0fpi_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpi_st_sh[ens_355]), err.(t0fpi_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpi_st_sh[ens_370]), err.(t0fpi_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpi_st_sh[ens_385]), err.(t0fpi_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot_SU2_pi(x_plot,uprm_st) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_pi(x_plot,uprm_st) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpi_ph_vec[2][1]), err(t0fpi_ph_vec[2][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.24, 0.33])
        setp(ax.get_yticklabels(),visible=false)

        subplot(133) # Create the 3rd axis of a 3x1 array of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpi_st_sh[ens_340]), err.(t0fpi_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpi_st_sh[ens_346]), err.(t0fpi_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpi_st_sh[ens_355]), err.(t0fpi_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpi_st_sh[ens_370]), err.(t0fpi_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpi_st_sh[ens_385]), err.(t0fpi_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        errorbar(value.(phi2_sh[ens_340]), value.(t0fpi_sh[ens_340]), err.(t0fpi_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fpi_sh[ens_346]), err.(t0fpi_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fpi_sh[ens_355]), err.(t0fpi_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fpi_sh[ens_370]), err.(t0fpi_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fpi_sh[ens_385]), err.(t0fpi_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        i = 1
        for ind in ind_sym
            list = [1,2,3,5]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_pi(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot_SU2_pi(x_plot,uprm_combined) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            list = [1,2,3,4]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_pi(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fpi_ph_vec[3][1]), err(t0fpi_ph_vec[3][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.24, 0.33])
        setp(ax.get_yticklabels(),visible=false)
        subplots_adjust(wspace=0.02)

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/SU2_fpi.pdf")
    ##

    #t0fk SU2
        uprm_combined = uprm_plot_SU2[3]
        uprm_st = uprm_plot_SU2[2]
        uprm = uprm_plot_SU2[1]

        fig = figure("SU2 kaon")
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        color_beta = ["rebeccapurple", "green", "blue", "darkorange", "red"]

        subplot(131) # Create the 1st axis of a 3x1 array of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fk_sh[ens_340]), err.(t0fk_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fk_sh[ens_346]), err.(t0fk_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fk_sh[ens_355]), err.(t0fk_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fk_sh[ens_370]), err.(t0fk_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fk_sh[ens_385]), err.(t0fk_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        ylabel(L"$\sqrt{8t_0}f_{K}$")
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot_SU2_k(x_plot,uprm) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_k(x_plot,uprm) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            #fill_between(x_plot[:,2], v-e, v+e, color_beta[i], alpha=0.75)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fk_ph_vec[1][1]), err(t0fk_ph_vec[1][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.3, 0.34])
        #legend()
        tight_layout()
            
        subplot(132) # Create the 2nd axis of a 3x1 arrax of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fk_st_sh[ens_340]), err.(t0fk_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fk_st_sh[ens_346]), err.(t0fk_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fk_st_sh[ens_355]), err.(t0fk_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fk_st_sh[ens_370]), err.(t0fk_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fk_st_sh[ens_385]), err.(t0fk_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot_SU2_k(x_plot,uprm_st) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_k(x_plot,uprm_st) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fk_ph_vec[2][1]), err(t0fk_ph_vec[2][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.3, 0.34])
        setp(ax.get_yticklabels(),visible=false)

        subplot(133) # Create the 3rd axis of a 3x1 array of axes
        errorbar(value.(phi2_sh[ens_340]), value.(t0fk_st_sh[ens_340]), err.(t0fk_st_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", mfc="none", label=L"$\beta=3.40$", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fk_st_sh[ens_346]), err.(t0fk_st_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", mfc="none", label=L"$\beta=3.46$", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fk_st_sh[ens_355]), err.(t0fk_st_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", mfc="none", label=L"$\beta=3.55$", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fk_st_sh[ens_370]), err.(t0fk_st_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", mfc="none", label=L"$\beta=3.70$", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fk_st_sh[ens_385]), err.(t0fk_st_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", mfc="none", label=L"$\beta=3.85$", color="red")
        errorbar(value.(phi2_sh[ens_340]), value.(t0fk_sh[ens_340]), err.(t0fk_sh[ens_340]), err.(phi2_sh[ens_340]), fmt="s", color="rebeccapurple")
        errorbar(value.(phi2_sh[ens_346]), value.(t0fk_sh[ens_346]), err.(t0fk_sh[ens_346]), err.(phi2_sh[ens_346]), fmt="o", color="green")
        errorbar(value.(phi2_sh[ens_355]), value.(t0fk_sh[ens_355]), err.(t0fk_sh[ens_355]), err.(phi2_sh[ens_355]), fmt="<", color="blue")
        errorbar(value.(phi2_sh[ens_370]), value.(t0fk_sh[ens_370]), err.(t0fk_sh[ens_370]), err.(phi2_sh[ens_370]), fmt=">", color="darkorange")
        errorbar(value.(phi2_sh[ens_385]), value.(t0fk_sh[ens_385]), err.(t0fk_sh[ens_385]), err.(phi2_sh[ens_385]), fmt="^", label=L"$\beta=3.85$", color="red")
        xlabel(L"$\phi_2$")
        i = 1
        for ind in ind_sym
            list = [1,2,3,4,7,8]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_k(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
        aux = model_plot_SU2_k(x_plot,uprm_combined) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        i = 1
        for ind in ind_sym
            list = [1,2,3,4,5,6]
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
            aux = model_plot_SU2_k(x_plot,uprm_combined[list]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        errorbar(value(phi2_ph), value(t0fk_ph_vec[3][1]), err(t0fk_ph_vec[3][1]), err(phi2_ph), fmt="x", label="ph. point", color="black")
        ax = gca()
        ax[:set_ylim]([0.3, 0.34])
        setp(ax.get_yticklabels(),visible=false)
        subplots_adjust(wspace=0.02)

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/SU2_fk.pdf")
    ##

    uprm_combined = uprm_plot[3]

    Fph = model_plot([0 ./ t0_sh ./ 8 [phi2_ph for i in 1:length(t0_sh)] phi4_sh], uprm_combined[[1,2,4]])
    Fphi2 = model_plot([0 ./ t0_sh ./ 8 phi2_sh phi4_sh], uprm_combined[[1,2,4]])
    y_aux = t0fpik_sh .- (Fphi2 .- Fph)

    Fph_st = model_plot([0 ./ t0_sh ./ 8 [phi2_ph for i in 1:length(t0_sh)] phi4_sh], uprm_combined[[1,2,3]])
    Fphi2_st = model_plot([0 ./ t0_sh ./ 8 phi2_sh phi4_sh], uprm_combined[[1,2,3]])
    y_aux_st = t0fpik_st_sh .- (Fphi2_st .- Fph_st)

    y_cont = t0fpik_ph_vec[3][1] 
    y_cont = t0fpik_ph[3] 
    delta_a = (y_aux[end-1] - y_cont) / err(y_cont); uwerr(delta_a); delta_a
    delta_a_st = (y_aux_st[end-1] - y_cont) / err(y_cont); uwerr(delta_a_st); delta_a_st

#==============================================================================================================================#

#============================== t0 sym à la Strassberger ======================================================================#

    y = sqrt.(t0_sh ./ t0_sh_sym)
    x = [1 ./ (8 .* t0_sh) phi2_sh phi4_sh phi2_sym]
    list = [2,3,5,7,8,9,10,12,13,14]
    y_aux = deepcopy(y)
    x_aux = deepcopy(x)
    y = y[list]
    x = x[list,:]
    Wm = inv(Symmetric(cov(y)))
    Wm = convert(Matrix{Float64}, Wm)

    function fun(x,p)
        return [sqrt(1 + p[1] * (x[i,2] - x[i,4])) for i in 1:length(x[:,1])]
    end
    uprm, chi_exp, chi2, pval_aux, doff = fit_alg(fun, value.(x), y, 1, diagm(diag(Wm)))

    ix = 3
    x_ph_aux = [0.0 8 * sqrt_t0_ph[ix] ^ 2 * Mpi ^ 2 / hc ^ 2 phi4_ph phi2_sym[1]]
    sqrt_t0_star = sqrt_t0_ph[ix] / fun(x_ph_aux,uprm)[1]; uwerr(sqrt_t0_star)
    details(sqrt_t0_star, string("syst chiral ",ix," 3rd")); sqrt(1-54/100) * err(sqrt_t0_star)

    R = 1 / fun(x_ph_aux,uprm)[1]

    a = sqrt_t0_ph[ix] * R ./ sqrt.(t0_sh_sym[ind_sym]); uwerr.(a)
    i = 1
    details(a[i], string("syst chiral ",ix," 3rd")); sqrt(1-26/100) * err(a[i])

    fig = figure()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 15
    x = deepcopy(x_aux); y = deepcopy(y_aux); uwerr.(y); uwerr.(x)
    errorbar(value.(x[:,2][ens_340]), value.(y[ens_340]), err.(y[ens_340]), err.(x[:,2][ens_340]), fmt="s", color="rebeccapurple")
    errorbar(value.(x[:,2][ens_346]), value.(y[ens_346]), err.(y[ens_346]), err.(x[:,2][ens_346]), fmt="o", color="green")
    errorbar(value.(x[:,2][ens_355]), value.(y[ens_355]), err.(y[ens_355]), err.(x[:,2][ens_355]), fmt="<", color="blue")
    errorbar(value.(x[:,2][ens_370]), value.(y[ens_370]), err.(y[ens_370]), err.(x[:,2][ens_370]), fmt=">", color="darkorange")
    errorbar(value.(x[:,2][ens_385]), value.(y[ens_385]), err.(y[ens_385]), err.(x[:,2][ens_385]), fmt="^", color="red")
    xlabel(L"$\phi_2$")
    ylabel(L"$\sqrt{t_0}/\sqrt{t_0^{\rm sym}}$")
    x_prime = [i for i in 0.0:0.01:0.75]
    x_plot = [0 .* x_prime x_prime [(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)]]
    aux = fun(x_plot,uprm) ; uwerr.(aux)
    v = value.(aux)
    e = err.(aux)
    fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.5)
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_sym.pdf")

    fb = BDIO_open("/home/asaez/cls_ens/results/t0_sym_combined.bdio", "w")
    write_uwreal(sqrt_t0_star ^ 2, fb, 1)
    BDIO_close!(fb)

#==============================================================================================================================#
