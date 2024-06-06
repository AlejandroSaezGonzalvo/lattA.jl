#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot, LinearAlgebra
using ADerrors: err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/plot.jl");

ens = ["H101", "H102r001", "H102r002", "H105", "H105r005", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N302", "J303", "E300", "J500", "J501"]
ens_av = ["H101", "H102", "H105", "H400", "D450", "N202", "N203", "N200", "D200", "E250", "N302", "J303", "E300", "J500", "J501"]
ens_340 = findall(x -> x in ["H101", "H102", "H105"], ens_av)
ens_346 = findall(x -> x in ["H400", "D450"], ens_av)
ens_355 = findall(x -> x in ["N202", "N203", "N200", "D200", "E250"], ens_av)
ens_370 = findall(x -> x in ["N300", "N302", "J303", "E300"], ens_av)
ens_385 = findall(x -> x in ["J500", "J501"], ens_av)#, "J501"]

#=========== 2203.08676 ===============#

w0_Mz = [
    uwreal([3.694,0.009], 1),
    uwreal([3.759,0.011], 2),
    uwreal([3.779,0.014], 3),
    uwreal([4.933,0.006], 4),
    uwreal([7.020,0.029], 5),
    uwreal([6.982,0.015], 6),
    uwreal([7.040,0.014], 7),
    uwreal([7.099,0.007], 8),
    uwreal([7.176,0.006], 9),
    uwreal([11.821,0.048], 10),
    uwreal([11.784,0.049], 11),
    uwreal([12.101,0.031], 12),
    uwreal([12.163,0.015], 13)
]

t0_Mz = [
    uwreal([2.847,0.005], 1),
    uwreal([2.881,0.006], 2),
    uwreal([2.889,0.007], 3),
    uwreal([3.696,0.003], 4),
    uwreal([5.165,0.012], 5),
    uwreal([5.146,0.006], 6),
    uwreal([5.164,0.006], 7),
    uwreal([5.179,0.003], 8),
    uwreal([5.203,0.002], 9),
    uwreal([8.544,0.019], 10),
    uwreal([8.526,0.019], 11),
    uwreal([8.621,0.010], 12),
    uwreal([8.622,0.006], 13)
]

uwerr.(w0_Mz)
uwerr.(t0_Mz)

#=========== read bdio ================#

obs = [Array{uwreal,1}() for i in 1:length(ens)]
for i in 1:length(ens)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/new_plateaux_noexp/unshifted/", ens[i], "_obs_wil_un.bdio"), "r")
    BDIO_seek!(fb); push!(obs[i], read_uwreal(fb))
    while BDIO_seek!(fb, 2) == true push!(obs[i], read_uwreal(fb)) end
    BDIO_close!(fb)
end
mpi = [obs[i][2] for i in 1:length(obs)]
mk = [obs[i][3] for i in 1:length(obs)]

ind = ensemble_inv["H102r002"]
mpi[ind-1] = plat_av(mpi, [ind-1,ind]); deleteat!(mpi, ind)
mk[ind-1] = plat_av(mk, [ind-1,ind]); deleteat!(mk, ind)

ind = ensemble_inv["H105r005"] - 1
mpi[ind-1] = plat_av(mpi, [ind-1,ind]); deleteat!(mpi, ind)
mk[ind-1] = plat_av(mk, [ind-1,ind]); deleteat!(mk, ind)

w0 = Array{uwreal,1}()
t0 = Array{uwreal,1}()
for id in ens
    obs = Array{uwreal,1}()
    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/w0_t0_", id, "_obs_wil_un.bdio"), "r")
    BDIO_seek!(fb); push!(obs, read_uwreal(fb)); BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb))
    BDIO_close!(fb)
    push!(w0,obs[1])
    push!(t0,obs[2])
end

t0[2] = plat_av(t0, [2,3])
t0[4] = plat_av(t0, [4,5])
w0[2] = plat_av(w0, [2,3])
w0[4] = plat_av(w0, [4,5])
deleteat!(t0,3)
deleteat!(t0,4)
deleteat!(w0,3)
deleteat!(w0,4)

phi4 = 8 * t0 .* (mk .^ 2 .+ 0.5 * mpi .^ 2)

uwerr.(phi4)
uwerr.(w0)
uwerr.(t0)

#============ plot ===================#

y_Mz= sqrt.(t0_Mz) ./ sqrt.(w0_Mz)
uwerr.(y_Mz)
x_Mz = 1 ./ t0_Mz
uwerr.(x_Mz)

y = sqrt.(t0) ./ sqrt.(w0)
uwerr.(y)
x = [1 ./ t0 phi4]
uwerr.(x)
C = Symmetric(cov(y))
W = convert(Matrix{Float64}, inv(C))

function model(x,p)
    return [p[1] + p[2] * x[i,1] + p[3] * x[i,1] ^ 2 for i in 1:length(x[:,1])] #+ p[4] * x[i,2]
end

uprm, chi_exp, chi2, pval_aux, doff = fit_alg(model, value.(x), y, 3, W)

Fph = model([0 ./ t0 [phi4_ph for i in 1:length(t0)]], uprm)
Fphi4 = model([0 ./ t0 phi4], uprm)
y_proj = y .- (Fphi4 .- Fph)
uwerr.(y_proj)

errorbar(value.(x[ens_340]), value.(y[ens_340]), err.(y[ens_340]), err.(x[ens_340]), fmt="s", label=L"$\beta=3.40$")
errorbar(value.(x[ens_346]), value.(y[ens_346]), err.(y[ens_346]), err.(x[ens_346]), fmt="<", label=L"$\beta=3.46$")
errorbar(value.(x[ens_355]), value.(y[ens_355]), err.(y[ens_355]), err.(x[ens_355]), fmt=">", label=L"$\beta=3.55$")
errorbar(value.(x[ens_370]), value.(y[ens_370]), err.(y[ens_370]), err.(x[ens_370]), fmt="^", label=L"$\beta=3.70$")
errorbar(value.(x[ens_385]), value.(y[ens_385]), err.(y[ens_385]), err.(x[ens_385]), fmt="o", label=L"$\beta=3.85$")

errorbar(value.(x_Mz), value.(y_Mz), err.(y_Mz), err.(x_Mz), fmt="x")

#errorbar(value.(x[ens_340]), value.(y_proj[ens_340]), err.(y_proj[ens_340]), err.(x[ens_340]), fmt="s")
#errorbar(value.(x[ens_346]), value.(y_proj[ens_346]), err.(y_proj[ens_346]), err.(x[ens_346]), fmt="<")
#errorbar(value.(x[ens_355]), value.(y_proj[ens_355]), err.(y_proj[ens_355]), err.(x[ens_355]), fmt=">")
#errorbar(value.(x[ens_370]), value.(y_proj[ens_370]), err.(y_proj[ens_370]), err.(x[ens_370]), fmt="^")
#errorbar(value.(x[ens_385]), value.(y_proj[ens_385]), err.(y_proj[ens_385]), err.(x[ens_385]), fmt="o")

x_plot = [i for i in 0.0:0.01:0.35]
aux = model([x_plot [phi4_ph for i in 1:length(x_plot)]], uprm) ; uwerr.(aux)
v = value.(aux)
e = err.(aux)
#fill_between(x_plot, v-e, v+e, color="gray", alpha=0.5)
legend()
xlabel(L"$a^2/t_0$")
ylabel(L"$\sqrt{t_0}/w_0$")
tight_layout()
