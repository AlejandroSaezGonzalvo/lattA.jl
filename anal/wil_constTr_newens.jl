import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "E300"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

#======== read correlators ===========#

if id == "N302"
    pp_sym, ap_sym = read_ens_TSM(path, ens) 
else
    pp_sym, ap_sym = get_corr_TSM_multichunks(path, ens)
end

#======== compute observables ========#

#tm = [[10], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
#tM = [[11], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
tm = [[10], collect(10:10:div(ens.T,2)-5)]
tM = [[ens.T-10], collect(ens.T-10:-10:div(ens.T,2)+5)]

if ens.id == "J501"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=[[10], [50,52,59,68,79]], tM=[[ens.T-10], [85,92,100,102,111,126,144]])
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=[[10], [50,52,59,68,79]], tM=[[ens.T-10], [85,92,100,102,111,126,144]])
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=[[10], [12,16,25]], tM=[[ens.T-10], [150,160,175]])
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=[[10], [12,16,25]], tM=[[ens.T-10], [150,160,175]])
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[ens.T-10], 96 .+ [10,15,20,25,40]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=true, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[ens.T-10], 96 .+ [10,15,20,25,40]])
elseif ens.id == "E300"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=[[10], [40,45,50,60,75]], tM=[[ens.T-10], [110,120,124,128]])
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=[[10], [35,40,55,65,90]], tM=[[ens.T-10], [105,112,129,139,154]])
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=[[10], [12,16,25]], tM=[[ens.T-10], [150,160,175]])
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=[[10], [12,16,25]], tM=[[ens.T-10], [150,160,175]])
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[ens.T-10], 96 .+ [10,15,20,25,40]])
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=true, wpm=wpm, tm=[[10], 96 .- [10,15,20,25,40]], tM=[[ens.T-10], 96 .+ [10,15,20,25,40]])
elseif ens.id == "N302"
    mpi = get_m(pp_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
    mk = get_m(pp_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
    m12 = get_mpcac(pp_sym[1], ap_sym[1], ens, "pion_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
    m13 = get_mpcac(pp_sym[2], ap_sym[2], ens, "kaon_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
    fpi = get_f_wil(pp_sym[1], ap_sym[1], mpi[1], ens, "pion_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
    fk = get_f_wil(pp_sym[2], ap_sym[2], mk[1], ens, "kaon_wil", pl=true, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "E250"
    ##TODO
    aux = 1
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
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=[[10], [10,20,30,40,50]], tM=[[ens.T-10], [90,100,110]])
end

#======== save BDIO ===================#

obs = [t0, mpi, mk, m12, m13, fpi, fk]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)

#============ get md ==================#

if md_meas == true
    phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
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
        push!(obs_md, s2 + v2)
    end

    ## now compute only strange derivatives wrt phi4 to interpolate as in Ben Strassberger's Thesis
    phi2 = 8 * t0 * mpi ^ 2
    t0fpik = sqrt(8 * t0) * 2/3 * (fk + 0.5 * fpi)
    for a in [phi4; phi2; t0fpik]
        md_s = [[md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(a, dSdm, YW, WY)]
        md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
        v2 = s2 = 0
        for i in 1:length(md_v)
            v2 += md_v[i][2]
            s2 += md_s[i][2]
        end
        s2 += md_s[end][2]
        push!(obs_md, s2 + v2)
    end
    obs_md[end-1] = obs_md[end-1] / obs_md[end-2]; obs_md[end] = obs_md[end] / obs_md[end-2]

    fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_md_wil.bdio"), "w")
    for i in 1:length(obs_md) write_uwreal(obs_md[i], fb, i) end
    BDIO_close!(fb)
end







