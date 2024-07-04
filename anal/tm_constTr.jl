#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/path_csv.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "N200"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = true

#======== read correlators ===========#

if ens.id in ["H102r001", "H102r002", "H105", "H105r005", "N203", "N200"]
    pp_sym, ap_sym, ppw, apw = read_ens_csv(ens)
    pp_sym_wil, ap_sym_wil, corr, corr_val, corrw, dSdm, w = read_ens_wil(path, ens, legacy=true)
else
    pp_sym, ap_sym, corrw, dSdm, w = read_ens_tm(path, ens, legacy=true)
end

#======== compute observables ========#

if ens.id in ["H102r001", "H102r002"]
    tm = [[10], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
    tM = [[11], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
else
    tm = [[10], collect(10:10:div(ens.T,2)-5)]
    tM = [[11], collect(ens.T-10:-10:div(ens.T,2)+5)]
end

mpi = Array{uwreal,1}()
m12 = Array{uwreal,1}()
fpi = Array{uwreal,1}()
mk = Array{uwreal,1}()
m13 = Array{uwreal,1}()
fk = Array{uwreal,1}()
m34 = Array{uwreal,1}()
for i in 1:15:length(pp_sym)
    for j in 0:2
        println(i, " ", j)
        mpi_aux = get_m(pp_sym[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(mpi, mpi_aux[1])
        m12_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(m12, m12_aux[1])
        m34_aux = get_mpcac(pp_sym[i+j+12], ap_sym[i+j+12], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(m34, m34_aux[1])
        if ens.id == "N200"
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[1], [42]], tM=[[10], [56]])
        elseif ens.id == "H102r002"
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[1], [24]], tM=[[10], [39]])
        elseif ens.id == "H102r001"
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[1], [20]], tM=[[10], [65]])
        else
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM, pl=false)
        end
        push!(fpi, fpi_aux[1])
    end
end
for i in 4:15:length(pp_sym)
    for j in 0:8
        println(i, " ", j)
        mk_aux = get_m(pp_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
        push!(mk, mk_aux[1])
        m13_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
        push!(m13, m13_aux[1])
        if ens.id == "N200"
            fk_aux = get_f_tm(pp_sym[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=[[1],[45,50,55]], tM=[[10],[108,110,117]])
        else
            fk_aux = get_f_tm(pp_sym[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
        end
        push!(fk, fk_aux[1])
    end
end

#======== compute t0/a² ===============#

if ens.id in ["H102r002"]
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm)
else
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM)
end

#======== save BDIO ===================#

obs = [mpi, mk, m12, m13, fpi, fk]
obs_str = ["mpi", "mk", "m12", "m13", "fpi", "fk"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "w")
    for i in 1:length(obs[j]) write_uwreal(obs[j][i], fb, i) end
    BDIO_close!(fb)
end

#============ get md ==================#

phi2 = 8 * t0 * mpi[1] ^ 2
phi4 = 8 * t0 * (mk[1] ^ 2 + 0.5 * mpi[1] ^ 2)
obs = [t0, phi2, phi4, sqrt(8 * t0) * m12[1], sqrt(8 * t0) * m13[1], sqrt(8 * t0) * fpi[1], sqrt(8 * t0) * fk[1], sqrt(8 * t0) * 2/3 * (fk[1] + 0.5 * fpi[1])]
if md_meas == true
    obs_md = Array{uwreal,1}()
    for a in obs
        if ens.id in ["H102r001", "H102r002", "H105", "H105r005", "N203", "N200"]
            md_s = [[md_sea(a, dSdm, ppw[i].obs, w) for i in 1:length(ppw)]; [md_sea(a, dSdm, apw[i].obs, w) for i in 1:length(apw)]; md_sea(a, dSdm, YW, WY)]
        else
            md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
        end
        s1 = s2 = 0
        for i in 1:length(md_s)
            s1 += md_s[i][1]
            s2 += md_s[i][2]
        end
        push!(obs_md, s2)
    end

    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_md_tm.bdio"), "w")
    for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
    BDIO_close!(fb)
end





