#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, LsqFit, LinearAlgebra

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");

err = ADerrors.err

obs = [Array{uwreal,1}() for i in 1:length(ensemble_old)]
der = [Array{uwreal,1}() for i in 1:length(ensemble_old)]
der_sea = [Array{uwreal,1}() for i in 1:length(ensemble_old)]
for i in 1:length(ensemble_old)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ensemble_old[i], "_md_wil.bdio"), "r") 
    BDIO_seek!(fb); push!(der[i], read_uwreal(fb)) 
    while BDIO_seek!(fb, 2) == true push!(der[i], read_uwreal(fb)) end 
    BDIO_close!(fb)

    fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ensemble_old[i], "_md_tm.bdio"), "r") 
    BDIO_seek!(fb); push!(der_sea[i], read_uwreal(fb)) 
    while BDIO_seek!(fb, 2) == true push!(der_sea[i], read_uwreal(fb)) end 
    BDIO_close!(fb)

    fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ensemble_old[i], "_obs_wil_un.bdio"), "r")
    BDIO_seek!(fb); push!(obs[i], read_uwreal(fb))
    while BDIO_seek!(fb, 2) == true push!(obs[i], read_uwreal(fb)) end
    BDIO_close!(fb)
end
t0 = [obs[i][1] for i in 1:length(obs)]
mpi = [obs[i][2] for i in 1:length(obs)]
mk = [obs[i][3] for i in 1:length(obs)]
phi2 = 8 * t0 .* mpi .^ 2
phi4 = 8 * t0 .* (mk .^ 2 .+ 0.5 * mpi .^ 2)

der_phi4 = [der[i][1] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)]
der_t0 = [der[i][2] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)]
der_phi2 = [der[i][3] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)]
der_t0m12 = [der[i][4] * beta_ZA[ens_db[ensemble_old[i]][3]] / beta_ZP[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_t0m13 = [der[i][5] * beta_ZA[ens_db[ensemble_old[i]][3]] / beta_ZP[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_t0fpi = [der[i][6] * beta_ZA[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_t0fk = [der[i][7] * beta_ZA[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_t0fpik = [der[i][8] * beta_ZA[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)

der_sea_phi2 = [der_sea[i][2] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)]
der_sea_phi4 = [der_sea[i][3] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)]
der_sea_m12 = [der_sea[i][4] * beta_ZA[ens_db[ensemble_old[i]][3]] / beta_ZP[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_sea_m13 = [der_sea[i][5] * beta_ZA[ens_db[ensemble_old[i]][3]] / beta_ZP[ens_db[ensemble_old[i]][3]] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_sea_t0fpi = [der_sea[i][6] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_sea_t0fk = [der_sea[i][7] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)
der_sea_t0fpik = [der_sea[i][8] for i in 1:length(der)] ./ [der[i][1] for i in 1:length(der)] ./ sqrt(8)

uwerr.(der_t0fpik)
uwerr.(der_t0fpi)
uwerr.(der_t0fk)
uwerr.(der_t0m12)
uwerr.(der_t0m13)
uwerr.(der_phi2)
uwerr.(der_t0)
uwerr.(der_sea_t0fpik)
uwerr.(der_sea_t0fpi)
uwerr.(der_sea_t0fk)
uwerr.(der_sea_phi2)
uwerr.(der_sea_phi4)
uwerr.(der_sea_m12)

function fun(x,p) 
    return [p[1] + p[2] * x[i,1] + p[3] * x[i,2] for i in 1:length(x[:,1])]
end

function fun_t0(x,p) 
    return [p[1] * (1 + p[2] * x[i,1] + p[3] / p[1] ^ 2) for i in 1:length(x[:,1])]
end

function funn(x,p) 
    return [p[1] + p[2] * x[i,1] + p[3] * x[i,1] ^ 2 + (p[4] + p[5] * x[i,1]) * x[i,2] for i in 1:length(x[:,1])]
end

## t0fpik

x = value.([phi2 1 ./ t0])
y = der_t0fpik
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0fpik, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0fpik

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}f_{\pi K}^{\rm W}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0fpi

x = value.([phi2 1 ./ t0])
y = der_t0fpi
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0fpi, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0fpi

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}f_{\pi}^{\rm W}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0fk

x = value.([phi2 1 ./ t0])
y = der_t0fk
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0fk, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0fk

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}f_{K}^{\rm W}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## phi2

x = [value.([phi2[2:5] 1 ./ t0[2:5]]); value.([phi2[8:10] 1 ./ t0[8:10]]); value.([phi2[12] 1 / t0[12]])]
y = [der_phi2[2:5]; der_phi2[8:10]; der_phi2[12]]
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_phi2, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_phi2

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\phi_2^{\rm W}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(der_phi2[1:5]), err.(der_phi2[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(der_phi2[6]), err(der_phi2[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(der_phi2[7:10]), err.(der_phi2[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(der_phi2[11:end]), err.(der_phi2[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0/a^2

x = value.([phi2 1 ./ t0])
y = der_t0
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{dt_0^{\rm W}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0/a^2 alternative

function fun_t0(x,p)
    f = [p[3] * (1 + p[1] * x[i] + p[2] / p[3] ^ 2) for i in 1:5]
    g = [p[4] * (1 + p[1] * x[i] + p[2] / p[4] ^ 2) for i in 6:6]
    h = [p[5] * (1 + p[1] * x[i] + p[2] / p[5] ^ 2) for i in 7:10]
    k = [p[6] * (1 + p[1] * x[i] + p[2] / p[6] ^ 2) for i in 11:12]
    return [f;g;h;k]
end

function fun_t0_plot(x,p)
    return [p[3] * (1 + p[1] * x[i] + p[2] / p[3] ^ 2) for i in 1:length(x)]
end

x = value.(phi2)
y = der_t0
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0_alt, chi_exp, chi2, pval_aux, doff = fit_alg(fun_t0, x, y, 6, Wm)
uprm = uprm_t0_alt

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{dt_0^{\rm W}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [i for i in auxx]
aux = fun_t0_plot(x_plot, uprm[[1,2,3]]) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [i for i in auxx]
aux = fun_t0_plot(x_plot, uprm[[1,2,4]]) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [i for i in auxx]
aux = fun_t0_plot(x_plot, uprm[[1,2,5]]) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [i for i in auxx]
aux = fun_t0_plot(x_plot, uprm[[1,2,6]]) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0fpik tm

x = value.([phi2 1 ./ t0])
y = der_sea_t0fpik
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0fpik_sea, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0fpik_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}f_{\pi K}^{\rm tm}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0fpi tm

x = value.([phi2 1 ./ t0])
y = der_sea_t0fpi
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0fpi_sea, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0fpi_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}f_{\pi}^{\rm tm}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0fk tm

x = value.([phi2 1 ./ t0])
y = der_sea_t0fk
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0fk_sea, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0fk_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}f_{K}^{\rm tm}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## phi2 tm

x = value.([phi2 1 ./ t0])
y = der_sea_phi2
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_phi2_sea, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_phi2_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\phi_2^{\rm tm}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## phi4 tm

x = value.([phi2 1 ./ t0])
y = der_sea_phi4
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_phi4_sea, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_phi4_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\phi_4^{\rm tm}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0m12 tm

x = value.([phi2 1 ./ t0])
y = der_sea_m12
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_m12_sea, chi_exp, chi2, pval_aux, doff = fit_alg(funn, x, y, 5, Wm)
uprm = uprm_m12_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}m_{12}^{\rm tm,R}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()


## t0m13 tm

x = value.([phi2 1 ./ t0])
y = der_sea_m13
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_m13_sea, chi_exp, chi2, pval_aux, doff = fit_alg(funn, x, y, 5, Wm)
uprm = uprm_m13_sea

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}m_{12}^{\rm tm,R}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

#=

## t0m12

x = [value.([phi2[2:5] 1 ./ t0[2:5]]); value.([phi2[8:10] 1 ./ t0[8:10]]); value.([phi2[12] 1 / t0[12]])]
y = [der_t0m12[2:5]; der_t0m12[8:10]; der_t0m12[12]]
#up = fit_routine(fun,x,y,3)
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0m12, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0m12

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}m_{12}^{\rm W,R}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(der_t0m12[1:5]), err.(der_t0m12[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(der_t0m12[6]), err(der_t0m12[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(der_t0m12[7:10]), err.(der_t0m12[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(der_t0m12[11:end]), err.(der_t0m12[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0m13

x = [value.([phi2[2:5] 1 ./ t0[2:5]]); value.([phi2[8:10] 1 ./ t0[8:10]]); value.([phi2[12] 1 / t0[12]])]
y = [der_t0m13[2:5]; der_t0m13[8:10]; der_t0m13[12]]
#up = fit_routine(fun,x,y,3)
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0m13, chi_exp, chi2, pval_aux, doff = fit_alg(fun, x, y, 3, Wm)
uprm = uprm_t0m13

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}m_{13}^{\rm W,R}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = fun(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(der_t0m13[1:5]), err.(der_t0m13[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(der_t0m13[6]), err(der_t0m13[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(der_t0m13[7:10]), err.(der_t0m13[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(der_t0m13[11:end]), err.(der_t0m13[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()
=#

## t0m12

x = value.([phi2 1 ./ t0])
y = der_t0m12
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0m12, chi_exp, chi2, pval_aux, doff = fit_alg(funn, x, y, 5, Wm)
uprm = uprm_t0m12

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}m_{12}^{\rm W,R}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

## t0m13

x = value.([phi2 1 ./ t0])
y = der_t0m13
Wm = inv(Symmetric(cov(y))); Wm = convert(Matrix{Float64}, Wm)
uprm_t0m13, chi_exp, chi2, pval_aux, doff = fit_alg(funn, x, y, 5, Wm)
uprm = uprm_t0m13

fig = figure()
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
xlabel(L"$\phi_2$")
ylabel(L"\frac{d\sqrt{t_0}m_{13}^{\rm W,R}}{d\phi_4^{\rm (s)}}")
auxx = collect(0.00:0.005:0.8)
x_plot = [[i for i in auxx] [1 / value(t0[1]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="rebeccapurple", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[6]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="green", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[7]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="blue", alpha=0.5)
x_plot = [[i for i in auxx] [1 / value(t0[11]) for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="darkorange", alpha=0.5)
x_plot = [[i for i in auxx] [0.0 for i in 1:length(auxx)]]
aux = funn(x_plot, uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
fill_between(x_plot[:,1], v-e, v+e, color="gray", alpha=0.5)
errorbar(value.(phi2[1:5]), value.(y[1:5]), err.(y[1:5]), fmt="x", color="rebeccapurple", label=L"\beta=3.40")
errorbar(value(phi2[6]), value(y[6]), err(y[6]), fmt="x", color="green", label=L"\beta=3.46")
errorbar(value.(phi2[7:10]), value.(y[7:10]), err.(y[7:10]), fmt="x", color="blue", label=L"\beta=3.55")
errorbar(value.(phi2[11:end]), value.(y[11:end]), err.(y[11:end]), fmt="x", color="darkorange", label=L"\beta=3.70")
legend()
tight_layout()

close("all")

der = [uprm_t0fpik; uprm_phi2; uprm_t0; uprm_t0fpik_sea; uprm_phi2_sea; uprm_phi4_sea; uprm_m12_sea; uprm_t0fpi; uprm_t0fk; uprm_t0m12; uprm_t0m13; uprm_t0fpi_sea; uprm_t0fk_sea; uprm_t0_alt; uprm_m13_sea]
fb = BDIO_open(string("/home/asaez/cls_ens/results/derivatives/der_1q_newrw.bdio"), "w")
for i in 1:length(der)
    write_uwreal(der[i], fb, i)
end
BDIO_close!(fb)
