#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "H105"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = true

if ens.id == "H102r001"
    tm = [[1], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
    tM = [[11], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
else
    tm = [[1], collect(10:10:div(ens.T,2)-5)]
    tM = [[11], collect(ens.T-10:-10:div(ens.T,2)+5)]
end
#=
w0 = get_w0(path, ens, [40,60], rw=true, wpm=wpm, tm=tm, tM=tM, pl=false)
if ens.id in ["H102r001", "H102r002"]
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm)
else
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM)
end

obs = [w0, t0]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/w0_t0_", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)
=#
#======== read correlators ===========#

pp_sym, ap_sym, corr, corr_val, corrw, dSdm, w = read_ens_wil(path, ens, legacy=true)

#======== compute observables ========#

if ens.id == "J303"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [52]], tM=[[10], [151]])
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "H105"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [28]], tM=[[10], [71]])
elseif ens.id == "H105r005"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [26]], tM=[[10], [63]])
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
else
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
end
if ens.id == "H105"
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [35]], tM=[[10], [76]])
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [36]], tM=[[10], [57]])
elseif ens.id == "H105r005"
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [35]], tM=[[10], [76]])
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [14]], tM=[[10], [80]])
else
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
end
if ens.id == "J303"
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[10], [50,55,60,65,75]], tM=[[11], [100,110,125,130,140]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "N200"
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [20]], tM=[[80], [90]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "N203"
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [30,40,50,60]], tM=[[80], [95,100]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [70,80,85]], tM=[[80], [95,98,100]])
elseif ens.id == "D200"
    #fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [20]], tM=[[80], [50]])
    #fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [30]], tM=[[80], [92]])
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [60,65]], tM=[[80], [70,80]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [30,40,60]], tM=[[80], [80,92]])
elseif ens.id == "H102r001"
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=[[1], [28,36,40]], tM=[[10], [46,53,65]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=[[1], [28,36,40]], tM=[[10], [46,53,65]])
else
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=false, wpm=wpm, tm=tm, tM=tM)
end

mpi, mk, m12, m13, fpi, fk = mpi[1], mk[1], m12[1], m13[1], fpi[1], fk[1]

#ZA = beta_ZA[ens.beta]
#bAtil = 1 + 0.0472 * (6 / ens.beta)
#fpi_un, fk_un = deepcopy(fpi), deepcopy(fk)
#fpi = ZA * (1 + bAtil * m12) * fpi
#fk = ZA * (1 + bAtil * m13) * fk 

#m12_I = (1 + beta_bap[ens.beta] * m12) * m12
#m13_I = (1 + beta_bap[ens.beta] * m13) * m13

#======== compute t0/a² ===============#

if ens.id in ["H102r002"]
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm)
else
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM)
end

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

#============ get md & dm =============#

phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
a = phi4
md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
v1 = v2 = s1 = s2 = 0
for i in 1:length(md_v)
    v1 += md_v[i][1]
    v2 += md_v[i][2]
    s1 += md_s[i][1]
    s2 += md_s[i][2]
end
s1 += md_s[end][1]
s2 += md_s[end][2]
dm = (phi4_ph - phi4) / (s2 + v2)

fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/dm_", ens.id, "_phi4=", round(value(phi4_ph), digits=5), ".bdio"), "w")
write_uwreal(dm, fb, 1)
BDIO_close!(fb)

phi2 = 8 * t0 * mpi ^ 2
phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
obs = [t0, phi2, phi4, sqrt(8 * t0) * m12, sqrt(8 * t0) * m13, sqrt(8 * t0) * fpi, sqrt(8 * t0) * fk, sqrt(8 * t0) * 2/3 * (fk + 0.5 * fpi)]
if md_meas == true
    obs_md = Array{uwreal,1}()
    for a in [phi4; obs]
        md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
        md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
        v1 = v2 = s1 = s2 = 0
        for i in 1:length(md_v)
            v1 += md_v[i][1]
            v2 += md_v[i][2]
            s1 += md_s[i][1]
            s2 += md_s[i][2]
        end
        s1 += md_s[end][1]
        s2 += md_s[end][2]
        #push!(obs_md, 2*s1 + s2 + v1 + v2)
        push!(obs_md, s2 + v2)
    end

    fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_md_wil.bdio"), "w")
    for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
    BDIO_close!(fb)
end







