#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO
using lattA: get_m, EnsInfo, fve, read_ens_tm, get_mpcac

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/path_csv.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "D200"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = true

#======== read correlators ===========#

if ens.id in ["H102r001", "H102r002", "H105", "H105r005", "N203", "N200"]
    pp_sym, ap_sym, dSdm = read_ens_csv(ens)
else
    pp_sym, ap_sym, corrw, dSdm, w = read_ens_tm(path, ens, legacy=true)
end

#======== compute observables ========#

#tm = [[10], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
#tM = [[11], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
tm = [[10], collect(10:10:div(ens.T,2)-5)]
tM = [[11], collect(ens.T-10:-10:div(ens.T,2)+5)]

mpi = Array{uwreal,1}()
m12 = Array{uwreal,1}()
fpi = Array{uwreal,1}()
mk = Array{uwreal,1}()
m13 = Array{uwreal,1}()
fk = Array{uwreal,1}()
m34 = Array{uwreal,1}()
for i in 1:8:length(pp_sym)
    for j in 0:1
        println(i, " ", j)
        if ens.id == "D200"
            mpi_aux = get_m(pp_sym[i+j], ens, "pion_tm", wpm=wpm, tm=[[10], [20,30,36,40,60]], tM=[[11], [80,85,100]])
            push!(mpi, mpi_aux[1])
            m12_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "pion_tm", wpm=wpm, tm=[[1], [60,62,65]], tM=[[11], [70,75,80]])
            push!(m12, m12_aux[1])
            m34_aux = get_mpcac(pp_sym[i+j+6], ap_sym[i+j+6], ens, "ss_tm", wpm=wpm, tm=[[10], [30,40,50]], tM=[[11], [70,80,90,110]])
            push!(m34, m34_aux[1])
            #fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[8,10,20,25,30,35,40,45,50,60,65,70,75], [8,10,20,25,30,35,40,45,50,60,65,70,75]], tM=[[90,95,100,105,110,115,120], [90,95,100,105,110,115,120]])
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[1], [40,50,60]], tM=[[10], [70,80,90]])
            push!(fpi, fpi_aux[1])
        else
            mpi_aux = get_m(pp_sym[i+j], ens, "pion_tm", wpm=wpm, tm=[[10], [50,59,70,84]], tM=[[11], [125,135,137]])
            push!(mpi, mpi_aux[1])
            m12_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "pion_tm", wpm=wpm, tm=[[10], [20,30,40]], tM=[[11], [70,80,90,110]])
            push!(m12, m12_aux[1])
            m34_aux = get_mpcac(pp_sym[i+j+6], ap_sym[i+j+6], ens, "ss_tm", wpm=wpm, tm=[[10], [30,40,50]], tM=[[11], [70,80,90,110]])
            push!(m34, m34_aux[1])
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[1], [40,45,55,60]], tM=[[10], [140,145,150]])
            push!(fpi, fpi_aux[1])
        end
    end
end
for i in 3:8:length(pp_sym)
    for j in 0:3
        println(i, " ", j)
        if ens.id == "D200"
            mk_aux = get_m(pp_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=[[10], [20,30,33,40,60]], tM=[[11], [80,84,100]])
            push!(mk, mk_aux[1])
            m13_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=[[10], [45,50,55,60]], tM=[[11], [70,75,80]])
            push!(m13, m13_aux[1])
            #fk_aux = get_f_tm(pp_sym[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=[[40,45,50,55], [40,45,50,55]], tM=[[92,95,100,105], [92,95,100,105]])
            fk_aux = get_f_tm(pp_sym[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=[[1], [40,50,60,70]], tM=[[10], [75,80,90,100,110,115]])
            push!(fk, fk_aux[1])
        else
            mk_aux = get_m(pp_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=[[1], [50,55,]], tM=[[11], [111,135,155]])
            push!(mk, mk_aux[1])
            m13_aux = get_mpcac(pp_sym[i+j], ap_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=[[10], [65,70,75]], tM=[[11], [120,130,135,140]])
            push!(m13, m13_aux[1])
            fk_aux = get_f_tm(pp_sym[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=[[1], [66,73,79]], tM=[[10], [120,125,140,160]])
            push!(fk, fk_aux[1])
        end
    end
end

#======== compute t0/a² ===============#

if ens.id in ["H102r001", "H102r002"]
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
        md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
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





