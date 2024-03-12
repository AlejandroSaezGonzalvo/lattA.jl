#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "E250"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

if id == "N302"
    pp_sym, ap_sym = read_ens_TSM(path, ens) 
else
    pp_sym, ap_sym = get_corr_TSM_multichunks(path, ens)
end

#======== Wilson =====================#

tm = [[10], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
tM = [[11], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
#tm = [[10], collect(10:10:div(ens.T,2)-5)]
#tM = [[ens.T-10], collect(ens.T-10:-10:div(ens.T,2)+5)]

if ens.id == "J501"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], [50,52,59,68,79]], tM=[[11], [85,92,100,102,111,126,144]])
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[10], [50,52,59,68,79]], tM=[[11], [85,92,100,102,111,126,144]])
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], [12,16,25]], tM=[[11], [150,160,175]])
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[10], [12,16,25]], tM=[[11], [150,160,175]])
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[11], 96 .+ [10,15,20,25,40]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[11], 96 .+ [10,15,20,25,40]])
elseif ens.id == "E300"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], [40,45,50,60,75]], tM=[[11], [110,120,124,128]])
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[10], [35,40,55,65,90]], tM=[[11], [105,112,129,139,154]])
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], [12,16,25]], tM=[[11], [150,160,175]])
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[10], [12,16,25]], tM=[[11], [150,160,175]])
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[11], 96 .+ [10,15,20,25,40]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[11], 96 .+ [10,15,20,25,40]])
elseif ens.id == "N302"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "E250"
    ##TODO

    tm = [[10], collect(10:10:div(ens.T,2)-5)]
    tM = [[94], [94]]

    mpi = get_m_pbc(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[5,10,15,20,25,30,40,50,60], tM=[96], method="cosh")
    mk = get_m_pbc(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[5,10,15,20,25,30,40,50,60], tM=[96], method="cosh")
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    fpi = get_f_wil_pbc(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[5,10,15,20,25,30,40,50,60], tM=[96])
    fk = get_f_wil_pbc(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[5,10,15,20,25,30,40,50,60], tM=[96])
end

mpi, mk, m12, m13, fpi, fk = mpi[1], mk[1], m12[1], m13[1], fpi[1], fk[1]
mpi, fpi, fk = fve(mpi, mk, fpi, fk, ens)

ZA = beta_ZA[ens.beta]
bAtil = 1 + 0.0472 * (6 / ens.beta)
fpi = ZA * (1 + bAtil * m12) * fpi
fk = ZA * (1 + bAtil * m13) * fk 

m12_I = (1 + beta_bap[ens.beta] * m12) * m12
m13_I = (1 + beta_bap[ens.beta] * m13) * m13

#======== compute t0/a² ===============#

if ens.id in ["E300", "J501", "E250"]
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "N302"
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=[[10], [10,20,30,40,50]], tM=[[11], [90,100,110]])
end

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

#=========== Wtm ======================#

if ens.id in ["E300", "J501"]
    pp_tm = pp_sym[10:end]
    ap_tm = ap_sym[10:end]
elseif ens.id in ["N302"]
    pp_tm = pp_sym[4:end]
    ap_tm = ap_sym[4:end]
end

mpi = Array{uwreal,1}()
m12 = Array{uwreal,1}()
fpi = Array{uwreal,1}()
mk = Array{uwreal,1}()
m13 = Array{uwreal,1}()
fk = Array{uwreal,1}()
m34 = Array{uwreal,1}()
for i in 1:8:length(pp_tm)
    for j in 0:1
        println(i, " ", j)
        if ens.id == "N302"
            mpi_aux = get_m(pp_tm[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
            push!(mpi, mpi_aux[1])
            m12_aux = get_mpcac(pp_tm[i+j], ap_tm[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
            push!(m12, m12_aux[1])
            m34_aux = get_mpcac(pp_tm[i+j+6], ap_tm[i+j+6], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
            push!(m34, m34_aux[1])
            fpi_aux = get_f_tm(pp_tm[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
            push!(fpi, fpi_aux[1])
        else
            mpi_aux = get_m(pp_tm[i+j], ens, "pion_tm", wpm=wpm, tm=tm_mpi[ens.id], tM=tM_mpi[ens.id])
            push!(mpi, mpi_aux[1])
            m12_aux = get_mpcac(pp_tm[i+j], ap_tm[i+j], ens, "pion_tm", wpm=wpm, tm=tm_m12[ens.id], tM=tM_m12[ens.id])
            push!(m12, m12_aux[1])
            m34_aux = get_mpcac(pp_tm[i+j+6], ap_tm[i+j+6], ens, "pion_tm", wpm=wpm, tm=tm_m12[ens.id], tM=tM_m12[ens.id])
            push!(m34, m34_aux[1])
            fpi_aux = get_f_tm(pp_tm[i+j], mpi[end], ens, "pion_tm", wpm=wpm, tm=tm_fpi[ens.id], tM=tM_fpi[ens.id])
            push!(fpi, fpi_aux[1])
        end
    end
end
for i in 3:8:length(pp_tm)
    for j in 0:3
        println(i, " ", j)
        if ens.id == "N302"
            mk_aux = get_m(pp_tm[i+j], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
            push!(mk, mk_aux[1])
            m13_aux = get_mpcac(pp_tm[i+j], ap_tm[i+j], ens, "pion_tm", wpm=wpm, tm=tm, tM=tM)
            push!(m13, m13_aux[1])
            fk_aux = get_f_tm(pp_tm[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=tm, tM=tM)
            push!(fk, fk_aux[1])
        else
            mk_aux = get_m(pp_tm[i+j], ens, "kaon_tm", wpm=wpm, tm=tm_mk[ens.id], tM=tM_mk[ens.id])
            push!(mk, mk_aux[1])
            m13_aux = get_mpcac(pp_tm[i+j], ap_tm[i+j], ens, "pion_tm", wpm=wpm, tm=tm_m13[ens.id], tM=tM_m13[ens.id])
            push!(m13, m13_aux[1])
            fk_aux = get_f_tm(pp_tm[i+j], mk[end], ens, "kaon_tm", wpm=wpm, tm=tm_fk[ens.id], tM=tM_fk[ens.id])
            push!(fk, fk_aux[1])
        end
    end
end

#=
c=0
for i in 1:length(mpi)
    for j in i+c:i+c+1
        global a, b, fk[j] = fve(mpi[i], mk[j], fpi[i], fk[j], ens)
    end
    mpi[i] = a
    fpi[i] = b
    c+=1
end
=#

#======== save BDIO ===================#

obs = [mpi, mk, m12, m13, fpi, fk]
obs_str = ["mpi", "mk", "m12", "m13", "fpi", "fk"]
for j in 1:length(obs_str)
    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_", obs_str[j], "_tm_un.bdio"), "w")
    for i in 1:length(obs[j]) write_uwreal(obs[j][i], fb, i) end
    BDIO_close!(fb)
end











