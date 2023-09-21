import Pkg
Pkg.activate("/home/alejandrosaez/University/PhD/codes/lattA")

using lattA

include("/home/alejandrosaez/University/PhD/codes/lattA/src/const.jl");

id = "H101"
ens = EnsInfo(id, ens_db[id])

path = joinpath(path, ens.id)
rwf = read_ms1.(path, v=vrw[id])

pp, ppw, w = get_corr(path, ens, "G5", "G5", rw=rwf, info=true)
pp_sym = [corr_sym(pp[i],pp[i+1],+1) for i in 1:2:length(pp)-1]

mpi = get_m()