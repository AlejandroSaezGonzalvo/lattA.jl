import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "J500"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

if id == "N302"
    pp_sym, ap_sym = read_ens_TSM(path, ens) 
else
    pp_sym, ap_sym = get_corr_TSM_multichunks(path, ens)
end

#=========== Wilson ==================#

#tm = [[10], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
#tM = [[11], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
tm = [[10], collect(10:10:div(ens.T,2)-5)]
tM = [[ens.T-10], collect(ens.T-10:-10:div(ens.T,2)+5)]

if ens.id == "J500"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, tm=tm, tM=tM)
    mk = mpi
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, tm=tm, tM=tM)
    m13 = m12
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[ens.T-10], 96 .+ [10,15,20,25,40]])
    fk = fpi
end

mpi, mk, m12, m13, fpi, fk = mpi[1], mk[1], m12[1], m13[1], fpi[1], fk[1]
mpi, fpi, fk = fve(mpi, mk, fpi, fk, ens)
fk = fpi ## need to "impose" this after fve in case of sym ens
mk = mpi

ZA = beta_ZA[ens.beta]
bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = ZA * (1 + bAtil * m12) * fpi
fk = ZA * (1 + bAtil * m13) * fk

m12_I = (1 + beta_bap[ens.beta] * m12) * m12

#======== compute t0/a² ===============#

t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM, pl=true)

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

#============== Wtm ===================#

pp_tm = pp_sym[4:end]
ap_tm = ap_sym[4:end]

mpi = Array{uwreal,1}()
m12 = Array{uwreal,1}()
fpi = Array{uwreal,1}()
mk = Array{uwreal,1}()
m13 = Array{uwreal,1}()
fk = Array{uwreal,1}()
m34 = Array{uwreal,1}()
for i in 1:length(pp_tm)
    println(i)
    mpi_aux = get_m(pp_tm[i], ens, "pion_tm", tm=tm, tM=tM, wpm=wpm)
    push!(mpi, mpi_aux[1])
    m12_aux = get_mpcac(pp_tm[i], ap_tm[i], ens, "pion_tm", tm=tm, tM=tM, wpm=wpm)
    push!(m12, m12_aux[1])
    fpi_aux = get_f_tm(pp_tm[i], mpi[i], ens, "pion_tm", tm=tm, tM=tM, wpm=wpm)
    push!(fpi, fpi_aux[1])
end
fk = fpi
mk = mpi

#=
for i in 1:length(mpi) mpi[i], fpi[i], fk[i] = fve(mpi[i], mk[i], fpi[i], fk[i], ens) end
fk = fpi ## need to "impose" this after fve in case of sym ens
mk = mpi
=#

#======== save BDIO ===================#

obs = [mpi, m12, fpi]
obs_str = ["mpi", "m12", "fpi"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "w")
    for i in 1:length(obs[j]) write_uwreal(obs[j][i], fb, i) end
    BDIO_close!(fb)
end