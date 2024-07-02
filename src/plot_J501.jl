id = "J501"
#import Pkg; Pkg.activate("/home/asaez/cls_ens/codes/lattA.jl")

using Revise, lattA, juobs, ADerrors, BDIO

include("/home/asaez/cls_ens/codes/lattA.jl/src/const.jl");
include("/home/asaez/cls_ens/codes/lattA.jl/src/in.jl");

#id_ind = parse(Int64, ARGS[1])
#id = ensemble[id_ind]
#id = "E300"
ens = EnsInfo(id, ens_db[id])

path = "/home/asaez/cls_ens/data"

const md_meas = false

tm = [[1], collect(div(ens.T,3)-4:div(ens.T,3)+4)]
tM = [[80], collect(div(2*ens.T,3)-4:div(2*ens.T,3)+4)]
#tm = [[1], collect(10:10:div(ens.T,2)-5)]
#tM = [[80], collect(ens.T-10:-10:div(ens.T,2)+5)]
if ens.id == "E250"
    tm = [[10], collect(10:10:div(ens.T,2)-5)]
    tM = [[94], [94]]
elseif ens.id == "D450"
    tm = [[10], collect(10:10:div(ens.T,2)-5)]
    tM = [[62], [62]]
end
#=
w0 = get_w0(path, ens, [40,60], rw=true, wpm=wpm, tm=tm, tM=tM, pl=false)
if ens.id in ["E300", "J501", "E250", "D450"]
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=tm, tM=tM)
elseif ens.id == "N302"
    t0, YW, WY = get_t0(path, ens, [40,60], rw=true, info=true, wpm=wpm, tm=[[10], [10,20,30,40,50]], tM=[[11], [90,100,110]])
end

obs = [w0, t0]
fb = BDIO_open(string("/home/asaez/cls_ens/results/unshifted/w0_t0_", ens.id, "_obs_wil_un.bdio"), "w")
for i in 1:length(obs) write_uwreal(obs[i], fb, i) end
BDIO_close!(fb)
=#
#======== read correlators ===========#

if id == "N302"
    pp_sym, ap_sym = read_ens_TSM(path, ens) 
else
    pp_sym, ap_sym = get_corr_TSM_multichunks(path, ens)
end

corr = pp_sym[1]

corr_d = corr.obs
m_dat = 0.5 .* log.((corr_d[2:end-2] ./ corr_d[3:end-1]) .^ 2)
y0 = corr.y0
T = length(corr_d) - 1 - y0

err = ADerrors.err

isnothing(wpm) ? uwerr(m) : uwerr(m, wpm)                       
	    isnothing(wpm) ? uwerr.(m_dat) : [uwerr(m_dat[i], wpm) for i in 1:length(m_dat)]
      	isnothing(wpm) ? uwerr.(m_i) : [uwerr(m_i[i], wpm) for i in 1:length(m_i)]
        v = value(m)
        e = err(m)

        fig = figure("eff")
        errorbar(1:length(m_dat), value.(m_dat), err.(m_dat), fmt="x", color="black")
        ylabel(L"$am_\mathrm{eff}$")
        xlabel(L"$x_0/a$")
        ylim(v-7*e, v+20*e)

        xlim(20,170)
        ylim(0.0625,0.08)

        v = 0.0661319
        e = 0.000154289

        fill_between(1:length(m_dat), v-e, v+e, color="green", alpha=0.5)

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/J501_meff.pdf")



corr = pp_sym[10]

corr_d = corr.obs
m_dat = 0.5 .* log.((corr_d[2:end-2] ./ corr_d[3:end-1]) .^ 2)
y0 = corr.y0
T = length(corr_d) - 1 - y0

err = ADerrors.err

isnothing(wpm) ? uwerr(m) : uwerr(m, wpm)                       
	    isnothing(wpm) ? uwerr.(m_dat) : [uwerr(m_dat[i], wpm) for i in 1:length(m_dat)]
      	isnothing(wpm) ? uwerr.(m_i) : [uwerr(m_i[i], wpm) for i in 1:length(m_i)]
        v = value(m)
        e = err(m)

        fig = figure("eff")
        errorbar(1:length(m_dat), value.(m_dat), err.(m_dat), fmt="x", color="black")
        ylabel(L"$am_\mathrm{eff}$")
        xlabel(L"$x_0/a$")
        ylim(v-7*e, v+20*e)

        xlim(20,170)
        ylim(0.0625,0.08)

        v = 0.06306933078942167
        e = 0.00013231014870663114

        fill_between(1:length(m_dat), v-e, v+e, color="green", alpha=0.5)

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/J501_meff_tm.pdf")



        corr_pp = pp_sym[1]
        corr_ap = ap_sym[1]
        impr = true

        ap_dat = -corr_ap.obs[2:end-1] 
    pp_dat = corr_pp.obs[2:end-1] 
    der_ap = (ap_dat[3:end] .- ap_dat[1:end-2]) / 2
    if impr == true
        ca = ens.ca
        der2_pp = pp_dat[1:end-2] + pp_dat[3:end] - 2*pp_dat[2:end-1]	
        der_ap = der_ap + ca*der2_pp
    end
    mpcac_dat = der_ap ./ (2*pp_dat[2:end-1])
    
    y0 = corr_pp.y0
    T = length(pp_dat) - 1 - y0

    isnothing(wpm) ? uwerr(mpcac) : uwerr(mpcac, wpm)                       
	    isnothing(wpm) ? uwerr.(mpcac_dat) : [uwerr(mpcac_dat[i], wpm) for i in 1:length(mpcac_dat)]
      	isnothing(wpm) ? uwerr.(mpcac_i) : [uwerr(mpcac_i[i], wpm) for i in 1:length(mpcac_i)]
        v = value(mpcac)
        e = err(mpcac)

        fig = figure()
        errorbar(1:length(mpcac_dat), value.(mpcac_dat), err.(mpcac_dat), fmt="x", color="black")
        ylabel(L"$am_\mathrm{pcac}$")
        xlabel(L"$x_0/a$")
        ylim(0.0026, 0.0031)
        xlim(1,188)

        v = 0.002740437602395615
        e = 3.6355305829617563e-6

        fill_between(1:length(mpcac_dat), v-e, v+e, color="green", alpha=0.5)

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/J501_mpcac.pdf")





        corr_pp = pp_sym[10]
        corr_ap = ap_sym[10]
        impr = true

        ap_dat = -corr_ap.obs[2:end-1] 
    pp_dat = corr_pp.obs[2:end-1] 
    der_ap = (ap_dat[3:end] .- ap_dat[1:end-2]) / 2
    if impr == true
        ca = ens.ca
        der2_pp = pp_dat[1:end-2] + pp_dat[3:end] - 2*pp_dat[2:end-1]	
        der_ap = der_ap + ca*der2_pp
    end
    mpcac_dat = der_ap ./ (2*pp_dat[2:end-1])
    
    y0 = corr_pp.y0
    T = length(pp_dat) - 1 - y0

    isnothing(wpm) ? uwerr(mpcac) : uwerr(mpcac, wpm)                       
	    isnothing(wpm) ? uwerr.(mpcac_dat) : [uwerr(mpcac_dat[i], wpm) for i in 1:length(mpcac_dat)]
      	isnothing(wpm) ? uwerr.(mpcac_i) : [uwerr(mpcac_i[i], wpm) for i in 1:length(mpcac_i)]
        v = value(mpcac)
        e = err(mpcac)

        fig = figure()
        errorbar(1:length(mpcac_dat), value.(mpcac_dat), err.(mpcac_dat), fmt="x", color="black")
        ylabel(L"$am_\mathrm{pcac}$")
        xlabel(L"$x_0/a$")
        ylim(0.00024, 0.0009)
        xlim(1,188)

        v = 0.00041515015836789
        e = 4.699215410379977e-6

        fill_between(1:length(mpcac_dat), v-e, v+e, color="green", alpha=0.5)

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/J501_mpcac_tm.pdf")

        tight_layout()






        corr_pp = pp_sym[1]
        corr_ap = ap_sym[1]
        impr = true

        ap_dat = -corr_ap.obs
    pp_dat = corr_pp.obs
    T = length(pp_dat)
    y0 = corr_pp.y0

    if impr == true
        ca = ens.ca
        der_pp = (pp_dat[3:end] .- pp_dat[1:end-2]) / 2
        ap_dat = ap_dat[2:end-1] + ca * der_pp
        aux = exp.((collect(1:T-2) .- (T-1)/2) .* [m for k = 1:T-2])
    else
        aux = exp.((collect(0:T-1) .- (T-1)/2) .* [m for k = 1:T])
    end
    R_dat = ap_dat .* aux ./ [((pp_dat[T-y0])^2)^(1/4) for k = 1:length(ap_dat)]
    #f_dat = [sqrt(2) * sqrt(R_dat[i] ^ 2) / sqrt(m) for i in 1:length(R_dat)]
    R_dat = sqrt.(R_dat .^ 2)

    isnothing(wpm) ? uwerr(R) : uwerr(R, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
      	isnothing(wpm) ? uwerr.(R_i) : [uwerr(R_i[i], wpm) for i in 1:length(R_i)]
        v = value(R)
        e = err(R)

        fig = figure("eff")
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_{\pi}$")
        xlabel(L"$x_0/a$")
        
        ylim(0.0062, 0.0072)
        xlim(10,170)

        v = 0.006457432033000142
        e = 4.728628484200814e-5

        fill_between(1:length(R_dat), v-e, v+e, color="green", alpha=0.5)

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/J501_R.pdf")



















        corr_pp = pp_sym[10]

        T = length(corr_pp.obs)
    y0 = corr_pp.y0
    mu = corr_pp.mu

    pp_dat = corr_pp.obs
    aux = exp.((collect(0:T-1) .- (T-1)/2 ) .* [m for k in 1:T])
    R_dat = pp_dat .* aux ./ [((pp_dat[T-y0])^2)^(1/4) for k = 1:T]


    isnothing(wpm) ? uwerr(R) : uwerr(R, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
      	isnothing(wpm) ? uwerr.(R_i) : [uwerr(R_i[i], wpm) for i in 1:length(R_i)]
        v = value(R)
        e = err(R)

        fig = figure("eff")
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_{\pi}$")
        xlabel(L"$x_0/a$")
        
        ylim(0.0062, 0.0072)
        xlim(10,170)

        v = 0.08024461099478472
        e = 0.0005714852320969807

        fill_between(1:length(R_dat), v-e, v+e, color="green", alpha=0.5)

        tight_layout()

        savefig("/home/asaez/cls_ens/codes/lattA.jl/plots/J501_R_tm.pdf")







        