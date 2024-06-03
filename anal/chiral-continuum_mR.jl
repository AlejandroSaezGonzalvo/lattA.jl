#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

#=
Wilson:
phi12 => E300, D450, ~H105 (down) are up, symmetric bad
phi13 => E300 (up), ~J303 (up), H105 (down), ~N200 (down), N302 (up), symmetric bad

MA:
phi12 => E250 (up), E300 (d), D200 (u), J303 (d), H105 (u), H101 (d), H400 (d), J501 (d)
phi13 => E300 (u), D459 (d), H105 (d), J303 (u), N302 (u), N202 (u), N300 (u), J500 (u), ~J501 (u)
=#

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, LsqFit, LinearAlgebra
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/chiral-continuum_fits.jl");

new = false
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

## JA
m12_sh_JA = [uwreal([0.008892,0.000113], "1"), uwreal([0.006440,0.000120], "2"), uwreal([0.004010,0.000149], "3"), uwreal([0.007815,0.000132], "4"), uwreal([0.006850,0.000087], "5"), uwreal([0.004796,0.000090], "6"), uwreal([0.003181,0.000083], "7"), uwreal([0.001625,0.000079], "8"), uwreal([0.005270,0.000085], "9"), uwreal([0.001919,0.000091], "10")]
m13_sh_JA = [uwreal([0.008892,0.000113], "1"), uwreal([0.010053,0.000121], "2"), uwreal([0.011210,0.000125], "3"), uwreal([0.007815,0.000132], "4"), uwreal([0.006850,0.000087], "5"), uwreal([0.007938,0.000088], "6"), uwreal([0.008649,0.000081], "7"), uwreal([0.009448,0.000083], "8"), uwreal([0.005270,0.000085], "9"), uwreal([0.007048,0.000101], "10")]
ZA = [[beta_ZA[3.4] for i in 1:3]; [beta_ZA[3.46] for i in 4:4]; [beta_ZA[3.55] for i in 5:8]; [beta_ZA[3.70] for i in 9:10]]
ZP = [[beta_ZP[3.4] for i in 1:3]; [beta_ZP[3.46] for i in 4:4]; [beta_ZP[3.55] for i in 5:8]; [beta_ZP[3.70] for i in 9:10]]
m12_sh_JA = m12_sh_JA .* ZA ./ ZP
m13_sh_JA = m13_sh_JA .* ZA ./ ZP

#============================== read BDIO =====================================================================================#

    obs = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    obs_tm_sh = [Array{uwreal,1}() for i in 1:length(ens)]
    for i in 1:length(ens)
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
    phi12_st_JA = sqrt.(8 * t0_sh) .* m12_sh_JA
    phi13_st_JA = sqrt.(8 * t0_sh) .* m13_sh_JA
    phi12 = sqrt.(8 * t0_sh) .* mul
    phi13 = sqrt.(8 * t0_sh) .* (mul .+ mus) ./ 2
    phiK_sh = phi4_sh .- 0.5 * phi2_sh

    y1_st = phi12_st ./ phi13_st; uwerr.(y1_st)
    y1_st_JA = phi12_st_JA ./ phi13_st_JA; uwerr.(y1_st_JA)
    y1 = phi12 ./ phi13; uwerr.(y1)
    y2_st = 2 * phi13_st ./ phiK_sh + phi12_st ./ phi2_sh; uwerr.(y2_st)
    y2_st_JA = 2 * phi13_st_JA ./ phiK_sh + phi12_st_JA ./ phi2_sh; uwerr.(y2_st_JA)
    y2 = 2 * phi13 ./ phiK_sh + phi12 ./ phi2_sh; uwerr.(y2)

    #t0fpik_sh = sqrt.(8 * t0_sh) .* (2/3) .* (1/2 * fpi_matched .+ fk_matched)
    #t0fpik_st_sh = sqrt.(8 * t0_sh) .* (2/3) .* (1/2 * fpi_sh .+ fk_sh) 

    K_sh = 1 ./ (16 * pi ^ 2 .* t0fpik_sh .^ 2)
    K_st_sh = 1 ./ (16 * pi ^ 2 .* t0fpik_st_sh .^ 2)

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
    
    errorbar(value.(phi2_sh[ens_340] .+ 0.02), value.(y2_st_JA[ens_340]), err.(y2_st_JA[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
    errorbar(value.(phi2_sh[ens_346] .+ 0.02), value.(y2_st_JA[ens_346]), err.(y2_st_JA[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
    errorbar(value.(phi2_sh[ens_355] .+ 0.02), value.(y2_st_JA[ens_355]), err.(y2_st_JA[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
    errorbar(value.(phi2_sh[ens_370] .+ 0.02), value.(y2_st_JA[ens_370]), err.(y2_st_JA[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
    errorbar(value.(phi2_sh[ens_385] .+ 0.02), value.(y2_st_JA[ens_385]), err.(y2_st_JA[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
    
    ax = gca()
    #ax[:set_ylim]([0.285, 0.325])
    #legend()
    tight_layout()
    savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/2phi13_div_phiK_plus_phi12_div_phi2_w.pdf")

#==============================================================================================================================#

#============================== chiral & continuum limits =====================================================================#

    ## prepare data 
        TIC = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        chi = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        W = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        pval = [Array{Float64,1}(), Array{Float64,1}(), Array{Float64,1}()]
        y1_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        y2_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        phi12_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        phi13_ph_vec = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_phi_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_phi12_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_phi13_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_yphi_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
        uprm_y_plot = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]

        x = [1 ./ (8 .* t0_sh) phi2_sh phi4_sh phi2_sym K_sh]
        x_st = [1 ./ (8 .* t0_sh) phi2_sh phi4_sh phi2_sym K_st_sh]
        x_355 = x[[ens_346; ens_355; ens_370; ens_385],:]
        x_nosym = x[ind_nosym,:]
        x_nosym355 = x[ind_nosym355,:]
        x_mL_40 = x[ind_mL_40, :]
        x_mL_41 = x[ind_mL_41, :]
        x_mL_42 = x[ind_mL_42, :]
        x_3552 = x[[ens_355; ens_370; ens_385],:]
        x_phi204 = x[ind_phi204,:]
        x_st_355 = x_st[[ens_346; ens_355; ens_370; ens_385],:]
        x_st_nosym = x_st[ind_nosym,:]
        x_st_nosym355 = x_st[ind_nosym355,:]
        x_st_mL_40 = x_st[ind_mL_40, :]
        x_st_mL_41 = x_st[ind_mL_41, :]
        x_st_mL_42 = x_st[ind_mL_42, :]
        x_st_3552 = x_st[[ens_355; ens_370; ens_385],:]
        x_st_phi204 = x_st[ind_phi204,:]
        x_ph = [0.0 (phi2_ph) (phi4_ph) (phi2_sym_ph) 1 / (16 * pi ^ 2 * 8 * t0_ph * fpik_exp ^ 2)]

        phi12_355 = phi12[[ens_346; ens_355; ens_370; ens_385]]
        phi12_3552 = phi12[[ens_355; ens_370; ens_385]]
        phi12_nosym = phi12[ind_nosym]
        phi12_nosym355 = phi12[ind_nosym355]
        phi12_phi204 = phi12[ind_phi204]
        phi12_mL_42 = phi12[ind_mL_42]

        phi12_st_355 = phi12_st[[ens_346; ens_355; ens_370; ens_385]]
        phi12_st_3552 = phi12_st[[ens_355; ens_370; ens_385]]
        phi12_st_nosym = phi12_st[ind_nosym]
        phi12_st_nosym355 = phi12_st[ind_nosym355]
        phi12_st_phi204 = phi12_st[ind_phi204]
        phi12_st_mL_42 = phi12_st[ind_mL_42]

        phi13_355 = phi13[[ens_346; ens_355; ens_370; ens_385]]
        phi13_3552 = phi13[[ens_355; ens_370; ens_385]]
        phi13_nosym = phi13[ind_nosym]
        phi13_nosym355 = phi13[ind_nosym355]
        phi13_phi204 = phi13[ind_phi204]
        phi13_mL_42 = phi13[ind_mL_42]

        phi13_st_355 = phi13_st[[ens_346; ens_355; ens_370; ens_385]]
        phi13_st_3552 = phi13_st[[ens_355; ens_370; ens_385]]
        phi13_st_nosym = phi13_st[ind_nosym]
        phi13_st_nosym355 = phi13_st[ind_nosym355]
        phi13_st_phi204 = phi13_st[ind_phi204]
        phi13_st_mL_42 = phi13_st[ind_mL_42]

        y1_355 = y1[[ens_346; ens_355; ens_370; ens_385]]
        y1_3552 = y1[[ens_355; ens_370; ens_385]]
        y1_nosym = y1[ind_nosym]
        y1_nosym355 = y1[ind_nosym355]
        y1_phi204 = y1[ind_phi204]
        y1_mL_42 = y1[ind_mL_42]

        y1_st_355 = y1_st[[ens_346; ens_355; ens_370; ens_385]]
        y1_st_3552 = y1_st[[ens_355; ens_370; ens_385]]
        y1_st_nosym = y1_st[ind_nosym]
        y1_st_nosym355 = y1_st[ind_nosym355]
        y1_st_phi204 = y1_st[ind_phi204]
        y1_st_mL_42 = y1_st[ind_mL_42]

        y2_355 = y2[[ens_346; ens_355; ens_370; ens_385]]
        y2_3552 = y2[[ens_355; ens_370; ens_385]]
        y2_nosym = y2[ind_nosym]
        y2_nosym355 = y2[ind_nosym355]
        y2_phi204 = y2[ind_phi204]
        y2_mL_42 = y2[ind_mL_42]

        y2_st_355 = y2_st[[ens_346; ens_355; ens_370; ens_385]]
        y2_st_3552 = y2_st[[ens_355; ens_370; ens_385]]
        y2_st_nosym = y2_st[ind_nosym]
        y2_st_nosym355 = y2_st[ind_nosym355]
        y2_st_phi204 = y2_st[ind_phi204]
        y2_st_mL_42 = y2_st[ind_mL_42]

        cuts_y1 = [y1, y1_355, y1_3552, y1_nosym, y1_nosym355, y1_phi204, y1_mL_42]
        cuts_y1_st = [y1_st, y1_st_355, y1_st_3552, y1_st_nosym, y1_st_nosym355, y1_st_phi204, y1_st_mL_42]

        cuts_y2 = [y2, y2_355, y2_3552, y2_nosym, y2_nosym355, y2_phi204, y2_mL_42]
        cuts_y2_st = [y2_st, y2_st_355, y2_st_3552, y2_st_nosym, y2_st_nosym355, y2_st_phi204, y2_st_mL_42]

        cuts_phi12 = [phi12, phi12_355, phi12_3552, phi12_nosym, phi12_nosym355, phi12_phi204, phi12_mL_42]
        cuts_phi12_st = [phi12_st, phi12_st_355, phi12_st_3552, phi12_st_nosym, phi12_st_nosym355, phi12_st_phi204, phi12_st_mL_42]

        cuts_phi13 = [phi13, phi13_355, phi13_3552, phi13_nosym, phi13_nosym355, phi13_phi204, phi13_mL_42]
        cuts_phi13_st = [phi13_st, phi13_st_355, phi13_st_3552, phi13_st_nosym, phi13_st_nosym355, phi13_st_phi204, phi13_st_mL_42]
        
        cuts_x = [x, x_355, x_3552, x_nosym, x_nosym355, x_phi204, x_mL_42]
        cuts_x_st = [x_st, x_st_355, x_st_3552, x_st_nosym, x_st_nosym355, x_st_phi204, x_st_mL_42]

        set_y1 = [cuts_y1, cuts_y1_st]
        set_y2 = [cuts_y2, cuts_y2_st]
        set_phi12 = [cuts_phi12, cuts_phi12_st]
        set_phi13 = [cuts_phi13, cuts_phi13_st]
        set_x = [cuts_x, cuts_x_st]

        #Wm = [inv.(Symmetric.(cov.(set_y[i]))) for i in 1:length(set_y)]
        #Wm = [convert.(Matrix{Float64}, Wm[i]) for i in 1:length(set_y)]
    ##

    models_y1 = [y1_model3, y1_model35]
    models_y2 = [y2_model3, y2_model35]
    models_x = [x_model3, x_model35]
    param = [5, 7]
    
    for k in 1:length(set_y1)
        println("k = ", k)
        for i in 1:length(models_y1)
            println("i = ", i)
            for j in 1:length(cuts_y1)
                println("j = ", j)
                yy1 = set_y1[k][j]
                yy2 = set_y2[k][j]
                x = set_x[k][j]
                indx = findall(value.(yy1) .!= 1.0)
                n = length(x[:,2]) + param[i]
                uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_y2[i], models_x[i]], [value.(x[indx,:]), value.(x), value.(x)], [yy1[indx], yy2, x[:,2]], n)
                #n = param[i]
                #W = [convert(Matrix{Float64}, Symmetric(inv(cov(yy1[indx])))), convert(Matrix{Float64}, Symmetric(inv(cov(yy2))))]
                #uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_y2[i]], [value.(x[indx,:]), value.(x)], [yy1[indx], yy2], W, n)
                push!(TIC[k], chi2 - 2 * chi_exp)
                push!(pval[k], pval_aux)
                push!(y1_ph_vec[k], models_y1[i]([x_ph;x], uprm)[1])
                push!(y2_ph_vec[k], models_y2[i]([x_ph;x], uprm)[1])
                if i == j == 1 uprm_y_plot[k] = uprm end
            end
        end
    end

    models_y1 = [y1_model3, y1_model35]
    models_y2 = [phi13_model1, phi13_model15]
    models_x = [x_model3, x_model35]
    param = [5, 7]
    
    for k in 1:length(set_y1)
        for i in 1:length(models_y1)
            for j in 1:length(cuts_y1)
                yy1 = set_y1[k][j]
                yy2 = set_phi13[k][j]
                x = set_x[k][j]
                indx = findall(value.(yy1) .!= 1.0)
                n = length(x[:,2]) + param[i]
                uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_y2[i], models_x[i]], [value.(x[indx,:]), value.(x), value.(x)], [yy1[indx], yy2, x[:,2]], n)
                #n = param[i]
                #W = [convert(Matrix{Float64}, Symmetric(inv(cov(yy1[indx])))), convert(Matrix{Float64}, Symmetric(inv(cov(yy2))))]
                #uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_y2[i]], [value.(x[indx,:]), value.(x)], [yy1[indx], yy2], W, n)
                push!(TIC[k], chi2 - 2 * chi_exp)
                push!(pval[k], pval_aux)
                push!(y1_ph_vec[k], models_y1[i]([x_ph;x], uprm)[1])
                push!(y2_ph_vec[k], models_y2[i]([x_ph;x], uprm)[1])
                if i == j == 1 uprm_yphi_plot[k] = uprm end
            end
        end
    end

    models_y1 = [phi12_model1, phi12_model15]
    models_y2 = [phi13_model1, phi13_model15]
    models_x = [x_model1, x_model15]
    param = [5, 7]
    
    for k in 1:length(set_y1)
        for i in 1:length(models_y1)
            for j in 1:length(cuts_y1)
                yy1 = set_phi12[k][j]
                yy2 = set_phi13[k][j]
                x = set_x[k][j]
                indx = findall(value.(yy1) .!= 1.0)
                n = length(x[:,2]) + param[i]
                uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_y2[i], models_x[i]], [value.(x[indx,:]), value.(x), value.(x)], [yy1[indx], yy2, x[:,2]], n)
                push!(TIC[k], chi2 - 2 * chi_exp)
                push!(pval[k], pval_aux)
                push!(phi12_ph_vec[k], models_y1[i]([x_ph;x], uprm)[1])
                push!(phi13_ph_vec[k], models_y2[i]([x_ph;x], uprm)[1])
                if i == j == 1 uprm_phi_plot[k] = uprm end
            end
        end
    end

    models_y1 = [phi12_model5]
    models_x = [x_model5]
    param = [5]
    
    for k in 1:length(set_y1)
        for i in 1:length(models_y1)
            for j in 1:length(cuts_y1)
                yy1 = set_phi12[k][j]
                x = set_x[k][j]
                n = length(x[:,2]) + param[i]
                uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_x[i]], [value.(x), value.(x)], [yy1, x[:,2]], n)
                push!(TIC[k], chi2 - 2 * chi_exp)
                push!(pval[k], pval_aux)
                push!(phi12_ph_vec[k], models_y1[i]([x_ph;x], uprm)[1])
                if i == j == 1 uprm_phi12_plot[k] = uprm end
            end
        end
    end

    models_y1 = [phi13_model6]
    models_x = [x_model6]
    param = [5]
    
    for k in 1:length(set_y1)
        for i in 1:length(models_y1)
            for j in 1:length(cuts_y1)
                yy1 = set_phi13[k][j]
                x = set_x[k][j]
                n = length(x[:,2]) + param[i]
                uprm, chi2, chi_exp, pval_aux = fit_alg([models_y1[i], models_x[i]], [value.(x), value.(x)], [yy1, x[:,2]], n)
                push!(TIC[k], chi2 - 2 * chi_exp)
                push!(pval[k], pval_aux)
                push!(phi13_ph_vec[k], models_y1[i]([x_ph;x], uprm)[1])
                if i == j == 1 uprm_phi13_plot[k] = uprm end
            end
        end
    end

    phi13_y_ph_vec = [[y2_ph_vec[k][i] / (y1_ph_vec[k][i] / phi2_ph + 2 / (phi4_ph - 0.5 * phi2_ph)) for i in 1:length(y1_ph_vec[k])] for k in 1:2]
    phi12_y_ph_vec = [[y1_ph_vec[k][i] * phi13_y_ph_vec[k][i] for i in 1:length(y1_ph_vec[k])] for k in 1:2]
    m12_y_ph_vec = [phi12_y_ph_vec[k] ./ sqrt(8 * t0_ph) * hc for k in 1:length(phi12_y_ph_vec)]
    m13_y_ph_vec = [phi13_y_ph_vec[k] ./ sqrt(8 * t0_ph) * hc for k in 1:length(phi12_y_ph_vec)]
    [uwerr.(m12_y_ph_vec[k]) for k in 1:length(m12_y_ph_vec)]
    [uwerr.(m13_y_ph_vec[k]) for k in 1:length(m13_y_ph_vec)]

    m12_ph_vec = [phi12_ph_vec[k] ./ sqrt(8 * t0_ph) * hc for k in 1:length(phi12_ph_vec)]
    m13_ph_vec = [phi13_ph_vec[k] ./ sqrt(8 * t0_ph) * hc for k in 1:length(phi12_ph_vec)]
    [uwerr.(m12_ph_vec[k]) for k in 1:length(m12_ph_vec)]
    [uwerr.(m13_ph_vec[k]) for k in 1:length(m13_ph_vec)]



   



    #TIC = [TIC[k] .- minimum.(TIC)[k] for k in 1:length(TIC)]
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

#================================== plots =====================================================================================#

    color_beta = ["rebeccapurple", "green", "blue", "darkorange", "red"]
    uwerr.(phi2_sh)
    
    ## Wilson ratios
        uprm_plot = uprm_y_plot

        fig = figure()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{12}/\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(y1_st[ens_340]), err.(y1_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(y1_st[ens_346]), err.(y1_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(y1_st[ens_355]), err.(y1_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(y1_st[ens_370]), err.(y1_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(y1_st[ens_385]), err.(y1_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = y1_model3(x_plot,uprm_plot[2]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = y1_model3(x_plot,uprm_plot[2]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        subplot(122)
        xlabel(L"$\phi_2$")
        ylabel(L"$2\phi_{13}/\phi_K+\phi_{12}/\phi_2$")
        errorbar(value.(phi2_sh[ens_340]), value.(y2_st[ens_340]), err.(y2_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(y2_st[ens_346]), err.(y2_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(y2_st[ens_355]), err.(y2_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(y2_st[ens_370]), err.(y2_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(y2_st[ens_385]), err.(y2_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = y2_model3(x_plot,uprm_plot[2]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = y2_model3(x_plot,uprm_plot[2]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/mR_ratios_w.pdf")
    ##

    ## Wtm ratios
        uprm_plot = uprm_y_plot

        fig = figure()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{12}/\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(y1[ens_340]), err.(y1[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(y1[ens_346]), err.(y1[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(y1[ens_355]), err.(y1[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(y1[ens_370]), err.(y1[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(y1[ens_385]), err.(y1[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = y1_model3(x_plot,uprm_plot[1]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = y1_model3(x_plot,uprm_plot[1]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        subplot(122)
        xlabel(L"$\phi_2$")
        ylabel(L"$2\phi_{13}/\phi_K+\phi_{12}/\phi_2$")
        errorbar(value.(phi2_sh[ens_340]), value.(y2[ens_340]), err.(y2[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(y2[ens_346]), err.(y2[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(y2[ens_355]), err.(y2[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(y2[ens_370]), err.(y2[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(y2[ens_385]), err.(y2[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = y2_model3(x_plot,uprm_plot[1]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = y2_model3(x_plot,uprm_plot[1]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/mR_ratios_tm.pdf")
    ##

    ## Wilson y_ratio
        uprm_plot = uprm_yphi_plot

        fig = figure()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{12}/\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(y1_st[ens_340]), err.(y1_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(y1_st[ens_346]), err.(y1_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(y1_st[ens_355]), err.(y1_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(y1_st[ens_370]), err.(y1_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(y1_st[ens_385]), err.(y1_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = y1_model3(x_plot,uprm_plot[2]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = y1_model3(x_plot,uprm_plot[2]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        subplot(122)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(phi13_st[ens_340]), err.(phi13_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(phi13_st[ens_346]), err.(phi13_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(phi13_st[ens_355]), err.(phi13_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(phi13_st[ens_370]), err.(phi13_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(phi13_st[ens_385]), err.(phi13_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = phi13_model1(x_plot,uprm_plot[2]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = phi13_model1(x_plot,uprm_plot[2]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/mR_yratio_w.pdf")
    ##

    ## Wtm y_ratio
        uprm_plot = uprm_yphi_plot

        fig = figure()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{12}/\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(y1[ens_340]), err.(y1[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(y1[ens_346]), err.(y1[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(y1[ens_355]), err.(y1[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(y1[ens_370]), err.(y1[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(y1[ens_385]), err.(y1[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = y1_model3(x_plot,uprm_plot[1]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = y1_model3(x_plot,uprm_plot[1]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        subplot(122)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(phi13[ens_340]), err.(phi13[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(phi13[ens_346]), err.(phi13[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(phi13[ens_355]), err.(phi13[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(phi13[ens_370]), err.(phi13[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(phi13[ens_385]), err.(phi13[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = phi13_model1(x_plot,uprm_plot[1]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = phi13_model1(x_plot,uprm_plot[1]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/mR_yratio_tm.pdf")
    ##

    ## Wilson phi X
        uprm_plot = uprm_phi_plot

        fig = figure()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{12}$")
        errorbar(value.(phi2_sh[ens_340]), value.(phi12_st[ens_340]), err.(phi12_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(phi12_st[ens_346]), err.(phi12_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(phi12_st[ens_355]), err.(phi12_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(phi12_st[ens_370]), err.(phi12_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(phi12_st[ens_385]), err.(phi12_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = phi12_model1(x_plot,uprm_plot[2]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = phi12_model1(x_plot,uprm_plot[2]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        subplot(122)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(phi13_st[ens_340]), err.(phi13_st[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(phi13_st[ens_346]), err.(phi13_st[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(phi13_st[ens_355]), err.(phi13_st[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(phi13_st[ens_370]), err.(phi13_st[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(phi13_st[ens_385]), err.(phi13_st[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = phi13_model1(x_plot,uprm_plot[2]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = phi13_model1(x_plot,uprm_plot[2]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/mR_Tay_w.pdf")
    ##

    ## Wtm phi X
        uprm_plot = uprm_phi_plot

        fig = figure()
        rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
        rcParams["font.size"] = 20

        subplot(121)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{12}$")
        errorbar(value.(phi2_sh[ens_340]), value.(phi12[ens_340]), err.(phi12[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(phi12[ens_346]), err.(phi12[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(phi12[ens_355]), err.(phi12[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(phi12[ens_370]), err.(phi12[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(phi12[ens_385]), err.(phi12[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = phi12_model1(x_plot,uprm_plot[1]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = phi12_model1(x_plot,uprm_plot[1]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        subplot(122)
        xlabel(L"$\phi_2$")
        ylabel(L"$\phi_{13}$")
        errorbar(value.(phi2_sh[ens_340]), value.(phi13[ens_340]), err.(phi13[ens_340]), err.(phi2_sh[ens_340]), label="", fmt="s", color="rebeccapurple", capsize=10.0)
        errorbar(value.(phi2_sh[ens_346]), value.(phi13[ens_346]), err.(phi13[ens_346]), err.(phi2_sh[ens_346]), label="", fmt="o", color="green", capsize=10.0)
        errorbar(value.(phi2_sh[ens_355]), value.(phi13[ens_355]), err.(phi13[ens_355]), err.(phi2_sh[ens_355]), label="", fmt="<", color="blue", capsize=10.0)
        errorbar(value.(phi2_sh[ens_370]), value.(phi13[ens_370]), err.(phi13[ens_370]), err.(phi2_sh[ens_370]), label="", fmt=">", color="darkorange", capsize=10.0)
        errorbar(value.(phi2_sh[ens_385]), value.(phi13[ens_385]), err.(phi13[ens_385]), err.(phi2_sh[ens_385]), fmt="^", color="red", label=L"\beta=3.85")
        i = 1
        x_prime = [i for i in 0.01:0.05:0.85]
        for ind in ind_sym
            x_plot = [[value(1 / (8 * t0_sh[ind])) for i in 1:length(x_prime)] x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
            aux = phi13_model1(x_plot,uprm_plot[1]) ; uwerr.(aux)
            v = value.(aux)
            e = err.(aux)
            #plot(x_plot[:,2], v, color=color_beta[i], alpha=0.6, linestyle="--")
            fill_between(x_plot[:,2], v-e, v+e, color=color_beta[i], alpha=0.3)
            i += 1
        end
        x_prime = [i for i in 0.01:0.05:0.85]
        x_plot = [0.0 * x_prime x_prime [value(phi4_ph) for i in 1:length(x_prime)] [value(phi2_sym_ph) for i in 1:length(x_prime)] [value(x_ph[1,5]) for i in 1:length(x_prime)]]
        aux = phi13_model1(x_plot,uprm_plot[1]) ; uwerr.(aux)
        v = value.(aux)
        e = err.(aux)
        fill_between(x_plot[:,2], v-e, v+e, color="gray", alpha=0.75)
        ax = gca()
        #ax[:set_ylim]([0.285, 0.325])
        #legend()

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/mR_Tay_tm.pdf")
    ##

#==============================================================================================================================#
