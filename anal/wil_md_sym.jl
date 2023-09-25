import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
id = "H101"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

#======== read md, correlators & obs ===========#

pp, ppw, w = get_corr_wil(path, ens, "G5", "G5", rw=true, info=true, legacy=true);
ap, apw, w = get_corr_wil(path, ens, "G5", "G0G5", rw=true, info=true, legacy=true);
YW, WY = get_YW(path, ens, [40,60], rw=true);

pp_d1 = get_corr_wil(path, ens, "G5_d1", "G5_d1", rw=true, legacy=true);
pp_d2 = get_corr_wil(path, ens, "G5_d2", "G5_d2", rw=true, legacy=true);
ap_d1 = get_corr_wil(path, ens, "G5_d1", "G0G5_d1", rw=true, legacy=true);
ap_d2 = get_corr_wil(path, ens, "G5_d2", "G0G5_d2", rw=true, legacy=true);
dSdm = get_dSdm(path, ens)

pp_val = [[pp_d1[i], pp_d2[i]] for i in 1:length(pp_d1)];
ap_val = [[ap_d1[i], ap_d2[i]] for i in 1:length(ap_d1)];
corr = [[pp[i] for i in 1:length(pp)]; [ap[i] for i in 1:length(ap)]];
corr_val = [[pp_val[i] for i in 1:length(pp)]; [ap_val[i] for i in 1:length(ap)]];
corrw = [[ppw[i] for i in 1:length(pp)]; [apw[i] for i in 1:length(ap)]];

obs = Array{uwreal,1}()
fb = BDIO_open(string("/home/asaez/cls_ens/results/", ens.id, "_obs_wil_un.bdio"), "r")
BDIO_seek!(fb); push!(obs, read_uwreal(fb))
for i in 2:7 BDIO_seek!(fb, 2); push!(obs, read_uwreal(fb)) end
t0, mpi, mk, m12, m13, fpi, fk = obs

#======== get dm ===============================#

phi4 = 8 * t0 * (mk ^ 2 + 0.5 * mpi ^ 2)
phi4_s = [[md_sea(phi4, dSdm, corrw[i], w) for i in 1:length(corrw)]; md_sea(phi4, dSdm, YW, WY)]
phi4_v = [md_val(phi4, corr[i], corr_val[i]) for i in 1:length(corr)]
phi4_v1 = phi4_v2 = phi4_s1 = phi4_s2 = 0
for i in 1:length(phi4_v) 
    phi4_v1 += phi4_v[i][1] 
    phi4_v2 += phi4_v[i][2] 
    phi4_s1 += phi4_v[i][1]
    phi4_s2 += phi4_v[i][2] 
end
phi4_d = 2*phi4_s1 + phi4_s2 + phi4_v1 + phi4_v2
dm = (phi4_ph - phi4) / phi4_d			

obs_sh = Array{uwreal,1}()
for a in obs
    md_s = [md_sea(a, dSdm, corrw[i], w) for i in 1:length(corrw)]
    md_v = [md_val(a, corr[i], corr_val[i]) for i in 1:length(corr)]
    v1 = v2 = s1 = s2 = 0
    for i in 1:length(md_v)
        v1 += md_v[i][1]
        v2 += md_v[i][2]
        v1 += md_s[i][1]
        s2 += md_s[i][2]
    end
    md = 2*s1 + s2 + v1 + v2
    push!(obs_sh, a + dm * md)
end

			
