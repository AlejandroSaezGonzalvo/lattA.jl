#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO, PyPlot
err = ADerrors.err

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "H101"
ens = EnsInfo(id, ens_db_Mainz[id])

path = "/home/asaez/cls_ens/data"

#tm = [[2], [28]]
#tM = [[1], [150]]
tm = [[2], [21]]
tM = [[1], [74]]

plt.ioff()

w0, t0 = get_w0t0(path, ens, [25,70], rw=true, wpm=wpm, tm=tm, tM=tM, pl=true, npol=2, w0_guess=3.5)

obs = [w0, t0]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/new_w0_t0_", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

close("all")