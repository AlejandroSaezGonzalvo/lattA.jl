#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using lattA, juobs, ADerrors

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");

id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=true)
pp_sym = [corr_sym(pp[i],pp[i+1],+1) for i in 1:2:length(pp)-1];

mpi = get_m(pp_sym[1], id)