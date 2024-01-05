import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#=========== read bdio ================#

obs = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "r")
BDIO_seek!(fb); push!(obs, read_uwreal(fb))
for i in 2:3 BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb)) end
BDIO_close!(fb)
t0, mpi_w, mk_w = obs
phi4_w = 8 * t0 * (mk_w ^ 2 + 0.5 * mpi_w ^ 2)
phi2_w = 8 * t0 * mpi_w ^ 2

obs = [Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}(), Array{uwreal,1}()]
obs_str = ["phi2", "phi4", "m12", "fpik"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_", obs_str[j], "_tm_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "r")
    BDIO_seek!(fb); push!(obs[j], read_uwreal(fb))
    while BDIO_seek!(fb, 2) == true push!(obs[j], read_uwreal(fb)) end 
    BDIO_close!(fb)
end
phi2 = obs[1]
phi4 = obs[2]
m12 = obs[3]
fpik = obs[4]

#========= match & full twist =========#

y = [m12; phi4]
target_m12 = 0 .* m12
target = [target_m12; [phi4_ph for i in 1:length(phi4)]]
y = y .- target
uwerr.(y)
W = 1 ./ err.(y) .^2

kappa, mul = ens_kappa[id], ens_mul[id]
x = [[[kappa[1] for i in 1:3]; [kappa[2] for i in 1:3]; [kappa[3] for i in 1:3]] [mul; mul; mul]]

up, chi2, chi_exp, pv = fit_alg(match_sym,x,y,6,[kappa[2], mul[2]],wpm=wpm) ##kappa->up[1], mul->up[2]