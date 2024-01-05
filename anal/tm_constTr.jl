import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "D200"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

pp, ppw, w = get_corr_tm(path, ens, "G5", "G5", rw=true, info=true, legacy=true);
pp_sym = [corr_sym(pp[i], pp[i+24], +1) for i in 1:24];
ap, apw, w = get_corr_tm(path, ens, "G5", "G0G5", rw=true, info=true, legacy=true);
ap_sym = [corr_sym(ap[i], ap[i+24], -1) for i in 1:24];

dSdm = get_dSdm(path, ens)

corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

#======== compute observables ========#

tm = [collect(10:5:div(ens.T,2)-5), collect(10:5:div(ens.T,2)-5)]
tM = [collect(div(ens.T,2)+5:5:ens.T-10), collect(div(ens.T,2)+5:5:ens.T-10)]

mpi = Array{uwreal,1}()
m12 = Array{uwreal,1}()
fpi = Array{uwreal,1}()
mk = Array{uwreal,1}()
m13 = Array{uwreal,1}()
fk = Array{uwreal,1}()
m34 = Array{uwreal,1}()
for i in 1:8:length(pp_sym)
    for j in 0:1
        mpi_aux = get_m(pp_sym[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(mpi, mpi_aux[1])
        m12_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(m12, m12_aux[1])
        m34_aux = get_mpcac(pp_sym[i+j+6], ap_sym[i+j+6], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(m34, m34_aux[1])
        fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(fpi, fpi_aux[1])
    end
end
for i in 3:8:length(pp_sym)
    for j in 0:3
        mk_aux = get_m(pp_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
        push!(mk, mk_aux[1])
        m13_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(m13, m13_aux[1])
        fk_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
        push!(fk, fk_aux[1])
    end
end

mpi, fpi, fk = fve.(mpi, mk, fpi, fk, ens)

#======== mass shift with fit =========#

par = Array{uwreal,1}()
fb = BDIO_open("/home/asaez/cls_ens/results/der_1q.bdio", "r")
BDIO_seek!(fb); push!(par, read_uwreal(fb))
while BDIO_seek!(fb, 2) == true push!(par, read_uwreal(fb)) end 
BDIO_close!(fb)

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_tm_un.bdio"), "w")
for a in obs write_uwreal(a, fb, i) end
BDIO_close!(fb)

#============ get md ==================#

if md_meas == true
    obs_md = [Array{uwreal,1}() for a in obs]
    for i in 1:length(obs)
        for a in obs[i]
            md_s = [md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]
            s1 = s2 = 0
            for i in 1:length(md_s)
                s1 += md_s[i][1]
                s2 += md_s[i][2]
            end
            push!(obs_md[i], 2*s1 + s2)
        end
    end


    for j in 1:length(obs_st)
        fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_", obs_st[j],"_md_tm.bdio"), "w")
        for i in 1:length(obs[j]) write_uwreal(obs_md[j][i], fb, i) end
        BDIO_close!(fb)
    end
end





