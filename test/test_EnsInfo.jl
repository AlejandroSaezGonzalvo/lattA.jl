import Pkg
Pkg.activate("/home/alejandrosaez/University/PhD/codes/lattA.jl")

using lattA

include("/home/alejandrosaez/University/PhD/codes/lattA.jl/src/const.jl")

id = "H101"
ens = EnsInfo(id, ens_db[id])