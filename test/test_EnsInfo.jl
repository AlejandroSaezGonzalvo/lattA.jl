import Pkg
Pkg.activate("/home/alejandrosaez/University/PhD/codes/lattA")

using lattA

include("/home/alejandrosaez/University/PhD/codes/lattA/src/const.jl")

id = "H101"
ens = EnsInfo(id, ens_db[id])