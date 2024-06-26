#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/path_csv.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H105"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

if ens.id in ["H102r001", "H102r002", "H105", "H105r005", "N203", "N200"]
    pp_sym, ap_sym = read_ens_csv(ens)
else
    pp_sym, ap_sym, corrw, dSdm = read_ens_tm(path, ens, legacy=true)
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
            fpi_aux = get_f_tm(pp_sym[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=[[1], [42,50,55]], tM=[[10], [78,80,90,100]])
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
        if ens.id in ["H105", "H105r005"]
            mk_aux = get_m(pp_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=[[1], [23]], tM=[[10], [69]])
        else
            mk_aux = get_m(pp_sym[i+j], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
        end
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

#======== save BDIO ===================#

obs = [mpi, mk, m12, m13, fpi, fk]
obs_str = ["mpi", "mk", "m12", "m13", "fpi", "fk"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "w")
    for i in 1:length(obs[j]) write_uwreal(obs[j][i], fb, i) end
    BDIO_close!(fb)
end

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





