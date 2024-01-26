using juobs, LaTeXStrings, PyPlot, ADerrors
include("const.jl");

function get_m(corr::juobs.Corr, ens::EnsInfo, PS::String; 
    tm::Union{Vector{Vector{Int64}}, Nothing}=nothing, tM::Union{Vector{Vector{Int64}}, Nothing}=nothing, 
    pl::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    corr_d = corr.obs
    m_dat = 0.5 .* log.((corr_d[2:end-2] ./ corr_d[3:end-1]) .^ 2)
    y0 = corr.y0
    T = length(corr_d) - 1 - y0

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1
    
    guess = value(m_dat[Int64(round(T / 2))])
    m, syst, m_i, weight, pval = model_av([fit_exp, fit_const], m_dat, guess, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
    if pl == true
        isnothing(wpm) ? uwerr(m) : uwerr(m, wpm)                       
	isnothing(wpm) ? uwerr.(m_dat) : [uwerr(m_dat[i], wpm) for i in 1:length(m_dat)]
      	isnothing(wpm) ? uwerr.(m_i) : [uwerr(m_i[i], wpm) for i in 1:length(m_i)]
        v = value(m)
        e = err(m)
        
        figure()
	subplot(411)
        fill_between(1:length(m_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(m_dat), value.(m_dat), err.(m_dat), fmt="x", color="black")
        ylabel(L"$m_\mathrm{eff}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+10*e)

	subplot(412)
	ylabel(L"$m_i$")
	fill_between(1:length(m_i), v-e, v+e, color="green", alpha=0.5)
	errorbar(1:length(m_i), value.(m_i), err.(m_i), fmt="x", color="black")

	subplot(413)
	ylabel(L"$W_i$")
	bar(1:length(m_i), weight, color="green")

	subplot(414)
	ylabel("p-value")
	bar(1:length(m_i), pval, color="green")

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/m_",ens.id,"_",PS,".pdf"))
        close("all")
    end

    return m, syst, m_i, weight, pval
end

function get_mpcac(corr_pp::juobs.Corr, corr_ap::juobs.Corr, ens::EnsInfo, PS::String; 
    tm::Union{Vector{Vector{Int64}}, Nothing}=nothing, tM::Union{Vector{Vector{Int64}}, Nothing}=nothing, 
    impr::Bool=true, pl::Bool=false, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

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

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1
    
    guess = value(mpcac_dat[Int64(round(T / 2))])
    mpcac, syst, mpcac_i, weight, pval = model_av([fit_exp, fit_const], mpcac_dat, guess, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
    
    if pl == true
        isnothing(wpm) ? uwerr(mpcac) : uwerr(mpcac, wpm)                       
	    isnothing(wpm) ? uwerr.(mpcac_dat) : [uwerr(mpcac_dat[i], wpm) for i in 1:length(mpcac_dat)]
      	isnothing(wpm) ? uwerr.(mpcac_i) : [uwerr(mpcac_i[i], wpm) for i in 1:length(mpcac_i)]
        v = value(mpcac)
        e = err(mpcac)
        
        figure()
	    subplot(411)
        fill_between(1:length(mpcac_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(mpcac_dat), value.(mpcac_dat), err.(mpcac_dat), fmt="x", color="black")
        ylabel(L"$m_\mathrm{pcac}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+10*e)

        subplot(412)
        ylabel(L"$mpcac_i$")
        fill_between(1:length(mpcac_i), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(mpcac_i), value.(mpcac_i), err.(mpcac_i), fmt="x", color="black")

        subplot(413)
        ylabel(L"$W_i$")
        bar(1:length(mpcac_i), weight, color="green")

        subplot(414)
        ylabel("p-value")
        bar(1:length(mpcac_i), pval, color="green")

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/mpcac_",ens.id,"_",PS,".pdf"))
        close("all")
    end

    return mpcac, syst, mpcac_i, weight, pval        
end


function get_f_wil(corr_pp::juobs.Corr, corr_ap::juobs.Corr, m::uwreal, ens::EnsInfo, PS::String;
    tm::Union{Vector{Vector{Int64}}, Nothing}=nothing, tM::Union{Vector{Vector{Int64}}, Nothing}=nothing, 
    impr::Bool=true, pl::Bool=false, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

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

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-1-y0-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1
    
    guess = value(R_dat[Int64(round(T / 2))])
    R, syst, R_i, weight, pval = model_av([fit_exp, fit_const], R_dat, guess, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
    f = sqrt(2) * sqrt(R^2) / sqrt(m)
    if pl == true
        isnothing(wpm) ? uwerr(R) : uwerr(R, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
      	isnothing(wpm) ? uwerr.(R_i) : [uwerr(R_i[i], wpm) for i in 1:length(R_i)]
        v = value(R)
        e = err(R)
        
        figure()
	    subplot(411)
        fill_between(1:length(R_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+10*e)

        subplot(412)
        ylabel(L"$R_i$")
        fill_between(1:length(R_i), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_i), value.(R_i), err.(R_i), fmt="x", color="black")

        subplot(413)
        ylabel(L"$W_i$")
        bar(1:length(R_i), weight, color="green")

        subplot(414)
        ylabel("p-value")
        bar(1:length(R_i), pval, color="green")

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,".pdf"))
        close("all")
    end
    
    return f, syst, R_i, weight, pval    
end

function get_f_wil(corr_ppL::juobs.Corr, corr_ppR::juobs.Corr, corr_apL::juobs.Corr, corr_apR::juobs.Corr, 
    m::uwreal, ens::EnsInfo, PS::String; tm::Union{Vector{Int64}, Nothing}=nothing, 
    tM::Union{Vector{Int64}, Nothing}=nothing, impr::Bool=true, pl::Bool=false, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    pp_dat = (corr_ppL.obs + corr_ppR.obs[end:-1:1]) / 2
    T = length(corr_ppL.obs)
    y0 = corr_ppL.y0

    if impr == true
        ca = ens.ca
        der_ppL = (corr_ppL.obs[3:end] - corr_ppL.obs[1:end-2]) / 2
        der_ppR = (corr_ppR.obs[3:end] - corr_ppR.obs[1:end-2]) / 2
        apL_dat = -corr_apL.obs[2:end-1] + ca * der_ppL
        apR_dat = -corr_apR.obs[2:end-1] + ca * der_ppR
    else
        apL_dat = -corr_apL.obs[2:end-1]
        apR_dat = -corr_apR.obs[2:end-1]
    end
    f1 = [pp_dat[T - y0] for k = 1:length(apL_dat)]
    
    R_dat = ((apL_dat .* apR_dat ./ f1).^2).^(1/4)
    #f_dat = [sqrt(2) * sqrt(R_dat[i] ^ 2) / sqrt(m) for i in 1:length(R_dat)]

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-1-y0-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1
    
    guess = value(R_dat[Int64(round(T / 2))])
    R, syst, R_i, weight, pval = model_av([fit_exp, fit_const], R_dat, guess, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
    f = sqrt(2) * sqrt(R^2) / sqrt(m)
    if pl == true
        isnothing(wpm) ? uwerr(R) : uwerr(R, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
      	isnothing(wpm) ? uwerr.(R_i) : [uwerr(R_i[i], wpm) for i in 1:length(R_i)]
        v = value(R)
        e = err(R)
        
        figure()
	    subplot(411)
        fill_between(1:length(R_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+10*e)

        subplot(412)
        ylabel(L"$R_i$")
        fill_between(1:length(R_i), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_i), value.(R_i), err.(R_i), fmt="x", color="black")

        subplot(413)
        ylabel(L"$W_i$")
        bar(1:length(R_i), weight, color="green")

        subplot(414)
        ylabel("p-value")
        bar(1:length(R_i), pval, color="green")

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,".pdf"))
        close("all")
    end
        
    return f, syst, R_i, weight, pval
end

function get_f_tm(corr_pp::juobs.Corr, m::uwreal, ens::EnsInfo, PS::String;
    tm::Union{Vector{Vector{Int64}}, Nothing}=nothing, tM::Union{Vector{Vector{Int64}}, Nothing}=nothing, 
    pl::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    T = length(corr_pp.obs)
    y0 = corr_pp.y0
    mu = corr_pp.mu

    pp_dat = corr_pp.obs
    aux = exp.((collect(0:T-1) .- (T-1)/2 ) .* [m for k in 1:T])
    R_dat = pp_dat .* aux ./ [((pp_dat[T-y0])^2)^(1/4) for k = 1:T]

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-1-y0-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1
    
    guess = value(R_dat[Int64(round(T / 2))])

    R, syst, R_i, weight, pval = model_av([fit_exp, fit_const], R_dat, guess, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
    f = sqrt(2) * (mu[1] + mu[2]) * R / m^1.5
    if pl == true
        isnothing(wpm) ? uwerr(R) : uwerr(R, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
      	isnothing(wpm) ? uwerr.(R_i) : [uwerr(R_i[i], wpm) for i in 1:length(R_i)]
        v = value(R)
        e = err(R)
        
        figure()
	    subplot(411)
        fill_between(1:length(R_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+10*e)

        subplot(412)
        ylabel(L"$R_i$")
        fill_between(1:length(R_i), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_i), value.(R_i), err.(R_i), fmt="x", color="black")

        subplot(413)
        ylabel(L"$W_i$")
        bar(1:length(R_i), weight, color="green")

        subplot(414)
        ylabel("p-value")
        bar(1:length(R_i), pval, color="green")

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,".pdf"))
        close("all")
    end

    return f, syst, R_i, weight, pval
end

function get_f_tm(corr_ppL::juobs.Corr, corr_ppR::juobs.Corr, m::uwreal, ens::EnsInfo, PS::String;
    tm::Union{Vector{Vector{Int64}}, Nothing}=nothing, tM::Union{Vector{Vector{Int64}}, Nothing}=nothing, 
    pl::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    T = length(corr_ppL.obs)
    y0 = corr_ppL.y0
    mu = corr_ppL.mu

    ppL_dat = corr_ppL.obs
    ppR_dat = corr_ppR.obs
    f1 = [ppL_dat[T - y0] for k = 1:T]
    R_dat = ((ppL_dat .* ppR_dat ./ f1).^2).^(1/4)

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-1-y0-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1
    
    guess = value(R_dat[Int64(round(T / 2))])
    R, syst, R_i, weight, pval = model_av([fit_exp, fit_const], R_dat, guess, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
    f = sqrt(2) * (mu[1] + mu[2]) * R / m^1.5
    if pl == true
        isnothing(wpm) ? uwerr(R) : uwerr(R, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
      	isnothing(wpm) ? uwerr.(R_i) : [uwerr(R_i[i], wpm) for i in 1:length(R_i)]
        v = value(R)
        e = err(R)
        
        figure()
	    subplot(411)
        fill_between(1:length(R_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+10*e)

        subplot(412)
        ylabel(L"$R_i$")
        fill_between(1:length(R_i), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(R_i), value.(R_i), err.(R_i), fmt="x", color="black")

        subplot(413)
        ylabel(L"$W_i$")
        bar(1:length(R_i), weight, color="green")

        subplot(414)
        ylabel("p-value")
        bar(1:length(R_i), pval, color="green")

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,".pdf"))
        close("all")
    end

    return f, syst, R_i, weight, pval
end

function get_t0(path::String, ens::EnsInfo, plat::Vector{Int64};
    tm::Union{Vector{Vector{Int64}}, Nothing}=nothing, tM::Union{Vector{Vector{Int64}}, Nothing}=nothing, pl::Bool=false, 
    rw=false, npol::Int64=2, ws::ADerrors.wspace=ADerrors.wsg, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, 
    info::Bool=false)

    path_ms = joinpath(path, ens.id, "gf")
    path_ms = filter(x->occursin(".dat", x), readdir(path_ms, join=true))
    Y = read_ms.(path_ms, id=ens.id, dtr=ens.dtr) 

    nr = length(Y)
    Ysl = getfield.(Y, :obs)
    t = getfield.(Y, :t)
    t = t[1]
    id = getfield.(Y, :id)
    id = id[1]
    replica = size.(Ysl, 1)

    L = ens.L
    #T = length(Y[:,1]) - y0
    y0 = 1 ## assumes this is the case, hardcoded, some ensembles will not fulfil !
    println("WARNING!: make sure t_src is 1 in this ensemble")

    #Truncation
    if id in keys(ADerrors.wsg.str2id)
        n_ws = findfirst(x-> x == ws.str2id[id], ws.map_nob)
        if !isnothing(n_ws)
            ivrep_ws = ws.fluc[n_ws].ivrep

            if length(replica) != length(ivrep_ws)
                error("Different number of replicas")
            end

            for k = 1:length(replica)
                if replica[k] > ivrep_ws[k]
                    println("Automatic truncation in Ysl ", ivrep_ws[k], " / ", replica[k], ". R = ", k)
                    Ysl[k] = Ysl[k][1:ivrep_ws[k], :, :]
                elseif replica[k] < ivrep_ws[k]
                    error("Automatic truncation failed. R = ", replica[k], "\nTry using truncate_data!")
                end
            end
            replica = size.(Ysl, 1)
        end
    end
    
    tmp = Ysl[1]
    [tmp = cat(tmp, Ysl[k], dims=1) for k = 2:nr]
    nt0 = juobs.t0_guess(t, tmp, plat, L)
    xmax = size(tmp, 2)
    T = xmax - 1 - y0

    dt0 = iseven(npol) ? Int64(npol / 2) : Int64((npol+1) / 2)
    Y_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)

    if rw
        path_rw = joinpath(path, ens.id, "rwf")
        path_rw = filter(x->occursin(".dat", x), readdir(path_rw, join=true))
        if ens.id == "D200"
            rwf_1 = read_ms1.([path_rw[1]], v=ens.vrw)
            rwf_2 = read_ms1.([path_rw[2]], v=ens.vrw)
            rwf = [hcat(rwf_2[1],rwf_1[1])]
        else
            rwf = read_ms1.(path_rw, v=ens.vrw)
        end

        Ysl_r, W = juobs.apply_rw(Ysl, rwf)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica)
        WY_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)
    end
    for i = 1:xmax
        k = 1
        for j = nt0-dt0:nt0+dt0
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica)
                Y_aux[i, k] = WY_aux[i, k] / W_obs
            end
            k = k + 1
        end
    end

    isnothing(tm) ? tm = [[y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40], [i for i in Int(round(T / 3)):Int(round(T / 3))+10]] : tm=tm
    isnothing(tM) ? tM = [[T-10,T-15,T-20,T-25,T-30,T-35,T-40], [i for i in Int(round(2 * T / 3)):Int(round(2 * T / 3))+10]] : tM=tM
    @.fit_exp(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-x))
    @.fit_const(x,p) = p[1] * x ^ 0
    k1 = 5
    k2 = 1

    x = t[nt0-dt0:nt0+dt0]
    t2E_i = Array{uwreal,1}()
    syst_i = Array{uwreal,1}()
    for j in 1:2*dt0+1
        i = Int(round(dt0+0.5))
        dat = Y_aux[:,j].* x[j].^2 / L^3
        if j == i            
            t2E_aux, syst_aux, t2E_aux_i, weight, pval = model_av([fit_exp, fit_const], dat, 0.3, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
            if pl == true
                isnothing(wpm) ? uwerr(t2E_aux) : uwerr(t2E_aux, wpm)                       
                isnothing(wpm) ? uwerr.(dat) : [uwerr(dat[i], wpm) for i in 1:length(dat)]
                isnothing(wpm) ? uwerr.(t2E_aux_i) : [uwerr(t2E_aux_i[i], wpm) for i in 1:length(t2E_aux_i)]
                v = value(t2E_aux)
                e = err(t2E_aux)
                
                figure()
                subplot(411)
                fill_between(1:length(dat), v-e, v+e, color="green", alpha=0.5)
                errorbar(1:length(dat), value.(dat), err.(dat), fmt="x", color="black")
                ylabel(L"$t^2\left<E\right>$")
                xlabel(L"$x_0$")
                ylim(v-10*e, v+10*e)

                subplot(412)
                ylabel(L"$\left(t^2\left<E\right>\right)_i$")
                fill_between(1:length(t2E_aux_i), v-e, v+e, color="green", alpha=0.5)
                errorbar(1:length(t2E_aux_i), value.(t2E_aux_i), err.(t2E_aux_i), fmt="x", color="black")

                subplot(413)
                ylabel(L"$W_i$")
                bar(1:length(t2E_aux_i), weight, color="green")

                subplot(414)
                ylabel("p-value")
                bar(1:length(t2E_aux_i), pval, color="green")

                savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/t2E_",ens.id,".pdf"))
                close("all")
            end
        else 
            t2E_aux, syst_aux, t2E_aux_i, weight, pval = model_av([fit_exp, fit_const], dat, 0.3, tm=tm, tM=tM, k=[k1,k2], wpm=wpm)
        end
        push!(t2E_i, t2E_aux)
        push!(syst_i, syst_aux)
    end

    if t2E_i[end] < 0.30 || t2E_i[1] > 0.30
        println("WARNING!: extrapolating t0/a²")
    end

    #uwerr.(t2E_i)
    #uwerr.(syst_i)
    #println("t2E_i = ", t2E_i)
    #println("syst_i = ", syst_i)
  
    model(x, p) = get_model(x, p, npol)

    par, aux = fit_alg(model, x, t2E_i, npol)
    fmin(x, p) = model(x, p) .- 0.3
    t0 = root_error(fmin, t[nt0], par)

    if pl == true
        v = value.(t2E_i)
        e = err.(t2E_i)
        uwerr(t0)

        errorbar(x, v, e, fmt="x")
        errorbar(value(t0), 0.3, xerr=err(t0), fmt="x")
        ylabel(L"$t^2E$")
        xlabel(L"$t/a^2$")
        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/t0_",id,".pdf"))
        close("all")
    end

    if info && rw
        return t0, WY_aux, W_obs
    elseif info && !rw
        return t0, Y_aux
    else
        return t0
    end

end
