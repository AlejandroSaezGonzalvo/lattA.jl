#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#============== read md & obs ==================#

obs = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "r")
BDIO_seek!(fb); push!(obs, read_uwreal(fb))
for i in 2:7 BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb)) end
t0, mpi, mk, m12, m13, fpi, fk = obs

obs_md = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_md_wil.bdio"), "r")
BDIO_seek!(fb); push!(obs_md, read_uwreal(fb))
for i in 1:7 BDIO_seek!(fb, 2); push!(obs_md, read_uwreal(fb)) end
phi4_d = obs_md[1]; obs_md = obs_md[2:end]

#========== get dm & mass shift obs ============#

phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
dm = (phi4_ph - phi4) / phi4_d		
obs_sh = Array{uwreal,1}()
for i in 1:length(obs) push!(obs_sh, obs[i] + dm * obs_md[i]) end

#================= save BDIO ===================#

fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_sh_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
write_uwreal(dm, fb, 8)
BDIO_close!(fb)
			
