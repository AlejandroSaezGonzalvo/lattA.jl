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
        errorbar(1:length(m_dat), value.(m_dat), err.(m_dat), fmt="x", color="black")
        ylabel(L"$m_\mathrm{eff}$")
        xlabel(L"$x_0$")
        ylim(v-7*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/m_",ens.id,"_",PS,"_plat.pdf"))
        
        models = Array{String,1}()
        for i in 1:length(tm[1])
            for j in 1:length(tM[1])
            push!(models, string("[",tm[1][i],",",tM[1][j],"]"))
            end
        end
        for i in 1:length(tm[2])
            for j in 1:length(tM[2])
            push!(models, string("[",tm[2][i],",",tM[2][j],"]"))
            end
        end

        figure()
	    subplot(411)
        fill_between(1:length(m_dat), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(m_dat), value.(m_dat), err.(m_dat), fmt="x", color="black")
        ylabel(L"$m_\mathrm{eff}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+20*e)
        ax = gca()
        setp(ax.get_xticklabels(),visible=false)

        subplot(412)
        ylabel(L"$m_i$")
        fill_between(1:length(m_i), v-e, v+e, color="green", alpha=0.5)
        errorbar(1:length(m_i), value.(m_i), err.(m_i), fmt="x", color="black")
        ax = gca()
        ax[:set_xlim]([0, length(models)+1])
        setp(ax.get_xticklabels(),visible=false)

        subplot(413)
        ylabel(L"$W_i$")
        bar(1:length(m_i), weight, color="green")
        ax = gca()
        ax[:set_xlim]([0, length(models)+1])
        setp(ax.get_xticklabels(),visible=false)

        subplot(414)
        ylabel("p-value")
        bar(1:length(m_i), pval, color="green")
        ax = gca()
        ax[:set_xlim]([0, length(models)+1])
        plt.xticks(collect(1:length(m_i)), models)
        xticks(rotation=90)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/m_",ens.id,"_",PS,".pdf"))
        close("all")
    end

    return m, syst, m_i, weight, pval
end

function get_m_pbc(corr::juobs.Corr, ens::EnsInfo, PS::String; 
    tm::Union{Vector{Int64}, Nothing}=nothing, tM::Union{Vector{Int64}, Nothing}=nothing, 
    pl::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing,
    method::String="cosh")
    
    corr_d = corr.obs

    m_dat = 0.5 .* log.((corr_d[2:end-2] ./ corr_d[3:end-1]) .^ 2)
    if ens.id == "E250"
        global T = 192
        global Thalf = 96
    elseif ens.id == "D450"
        global T = 128
        global Thalf = 64
    end
    guess = value(m_dat[Int64(round(T / 3))])

    isnothing(tm) ? tm = [y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40] : tm=tm
    isnothing(tM) ? tM = [Thalf] : tM=tM

    aux = [corr_d[i] / corr_d[i+1] for i in 2:length(corr_d)-1]
    if ens.id == "E250"
        global T = 96
    elseif ens.id == "D450"
        global T = 64
    end
    @.fit_fun(x,p) = cosh(p[1] * (x-T)) / cosh(p[1] * (x+1-T))
    k1 = 1

    aux2 = corr_d
    global T2 = 192
    @.fit_fun2(x,p) = p[2] * exp(-p[1] * (x-1)) + p[2] * exp(-p[1] * (T2+1-x))
    k2 = 2

    if method == "cosh"
        m, syst, m_i, weight, pval = model_av(fit_fun, aux, guess, tm=tm, tM=tM, k=k1, wpm=wpm)
    elseif method == "corr"
        m, syst, m_i, weight, pval = model_av(fit_fun2, aux2, guess, tm=tm, tM=tM, k=k2, wpm=wpm)
    end

    if pl == true
        ##TODO
        bla = 1
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
        errorbar(1:length(mpcac_dat), value.(mpcac_dat), err.(mpcac_dat), fmt="x", color="black")
        ylabel(L"$m_\mathrm{pcac}$")
        xlabel(L"$x_0$")
        ylim(v-10*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/mpcac_",ens.id,"_",PS,"_plat.pdf"))
        
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
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-20*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.pdf"))
        
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
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-20*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.pdf"))
        
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

function get_f_wil_pbc(corr_pp::juobs.Corr, corr_ap::juobs.Corr, m::uwreal, ens::EnsInfo, PS::String;
    tm::Union{Vector{Int64}, Nothing}=nothing, tM::Union{Vector{Int64}, Nothing}=nothing, 
    impr::Bool=true, pl::Bool=false, method::String="ratio",
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    ap_dat = -corr_ap.obs
    pp_dat = corr_pp.obs
    if ens.id == "E250"
        global T = 192
        global Thalf = 96
    elseif ens.id == "D450"
        global T = 128
        global Thalf = 64
    end

    if impr == true
        ca = ens.ca
        der_pp = (pp_dat[3:end] .- pp_dat[1:end-2]) / 2
        ap_dat = ap_dat[2:end-1] + ca * der_pp
        pp_dat = pp_dat[2:end-1]
    end
    R_dat = ap_dat ./ sqrt.(pp_dat)

    isnothing(tm) ? tm = [y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40] : tm=tm
    isnothing(tM) ? tM = [Thalf] : tM=tM
    
    @.fit_fun(x,p) = p[1] * (-exp(-value(m) * (x-1)) + exp(-value(m) * (T+1-x))) / sqrt(exp(-value(m) * (x-1)) + exp(-value(m) * (T+1-x))) 
    @.fit_fun_ap(x,p) = p[1] * (-exp(-value(m) * (x-1)) + exp(-value(m) * (T+1-x)))
    @.fit_fun_pp(x,p) = p[1] * (exp(-value(m) * (x-1)) + exp(-value(m) * (T+1-x)))
    @.fit_const(x,p) = p[1] * x ^ 0
    k = 1

    function fit_corr_sim(x,p)
		f = [p[1] * (-exp(-value(m) * x[i]) + exp(-value(m) * (T - x[i]))) for i in 1:div(length(x),2)]
		g = [p[2] * (exp(-value(m) * x[i]) + exp(-value(m) * (T - x[i]))) for i in div(length(x),2)+1:length(x)]
		return [f;g]
	end
    
    if method == "ratio"
        R, syst, R_i, weight, pval = model_av(fit_fun, R_dat, .5, tm=tm, tM=tM, k=k, wpm=wpm)
        f_i = sqrt.(2 ./ [m for i in 1:length(R_i)]) .* R_i 
        f = sqrt(2 / m) * R 
    elseif method == "corr"
        cap, syst_ap, cap_i, weight_ap, pval_ap = model_av(fit_fun_ap, ap_dat, .5, tm=tm, tM=tM, k=k, wpm=wpm)
        cpp, syst_pp, cpp_i, weight_pp, pval_pp = model_av(fit_fun_pp, pp_dat, .5, tm=tm, tM=tM, k=k, wpm=wpm)
        f_i = sqrt.(2 ./ [m for i in 1:length(cap_i)]) .* cap_i ./ sqrt.(cpp_i) 
        f = sqrt(2 / m) * cap / sqrt(cpp)
        syst = [syst_ap, syst_pp]
        pval = [pval_ap, pval_pp]
        weight = [weight_ap, weight_pp]
    elseif method == "corr_sim"
        pval = Array{Float64,1}()
        TIC = Array{Float64,1}()
        cap_i = Array{uwreal,1}()
        cpp_i = Array{uwreal,1}()
        for i in tm
            y = [ap_dat[i:end]; pp_dat[i:end]]
            x = collect(i:length(ap_dat))
            x = [x; x]
            uwerr.(y)
            W = 1 ./ ADerrors.err.(y) .^ 2
            uprm, chi2, chi_exp, pv = fit_alg(fit_corr_sim, x, y, 2, [.5, .5])
            push!(cap_i, uprm[1])
            push!(cpp_i, uprm[2])
            push!(pval, pv)
            push!(TIC, chi2 - 2 * chi_exp)
        end
        weight = exp.(-0.5 * TIC) ./ sum(exp.(-0.5 * TIC))
        R_i = cap_i ./ sqrt.(cpp_i)
        R_av = sum(R_i .* weight)
        syst = sum(R_i .^ 2 .* weight) - R_av ^ 2
        f_i = sqrt.(2 ./ [m for i in 1:length(cap_i)]) .* R_i
        f = sqrt(2 / m) * R_av
    elseif method == "pcac"
        if impr == false
            error("pcac method implemented only for improved axial current, set impr=true")
        end
        cpp, syst_pp, cap_i, weight, pval = model_av(fit_fun_pp, pp_dat, .5, tm=tm, tM=tM, k=k, wpm=wpm)
        der_ap = (ap_dat[3:end] .- ap_dat[1:end-2]) / 2
        der2_pp = pp_dat[1:end-2] + pp_dat[3:end] - 2*pp_dat[2:end-1]	
        der_ap = der_ap + ca*der2_pp
        mpcac_dat = der_ap ./ (2*pp_dat[2:end-1])
        f_dat = [sqrt(8 / m ^ 3 * cpp) for i in 1:length(mpcac_dat)] .* mpcac_dat
        f, syst, f_i, weight, pval = model_av(fit_const, f_dat, .5, tm=tm, tM=tM, k=k, wpm=wpm)
    end

    if pl == true
        ##TODO
        bla = 1

        x = collect(1:length(R_dat))
        R_dat = R_dat ./ [(-exp(-value(m) * (x[i]-1)) + exp(-value(m) * (T+1-x[i]))) / sqrt(exp(-value(m) * (x[i]-1)) + exp(-value(m) * (T+1-x[i]))) for i in 1:length(x)] 

	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
        v = value(R_dat[40])
        e = err(R_dat[40])

        figure()
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-20*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.pdf"))
    end
    
    return f, syst, f_i, weight, pval    
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
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-20*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.pdf"))
        
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
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        ylim(v-20*e, v+20*e)
        tight_layout()

        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.pdf"))
        
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

function get_f_tm_pbc(corr_pp::juobs.Corr, m::uwreal, ens::EnsInfo, PS::String;
    tm::Union{Vector{Int64}, Nothing}=nothing, tM::Union{Vector{Int64}, Nothing}=nothing, 
    pl::Bool=false,
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)

    pp_dat = corr_pp.obs[2:end]
    mu = corr_pp.mu
    if ens.id == "E250"
        global T = 192
        global Thalf = 96
    elseif ens.id == "D450"
        global T = 128
        global Thalf = 64
    end

    isnothing(tm) ? tm = [y0+10,y0+15,y0+20,y0+25,y0+30,y0+35,y0+40] : tm=tm
    isnothing(tM) ? tM = [Thalf] : tM=tM
    
    @.fit_fun_pp(x,p) = p[1] * (exp(-value(m) * (x-1)) + exp(-value(m) * (T+1-x)))
    k = 1
    
    cpp, syst, cpp_i, weight, pval = model_av(fit_fun_pp, pp_dat, .5, tm=tm, tM=tM, k=k, wpm=wpm)
    f_i = [sqrt(2) * (mu[1] + mu[2]) / m ^ 1.5 for i in 1:length(cpp_i)] .* sqrt.(cpp_i) 
    f = sqrt(2) * (mu[1] + mu[2]) / m ^ 1.5 * sqrt(cpp)

    if pl == true
        ##TODO
        bla = 1

        x = collect(1:length(pp_dat))
        R_dat = pp_dat #./ [exp(-m * (x[i]-1)) + exp(-m * (T+1-x[i])) for i in 1:length(x)]
        vec = cpp .* [exp(-m * (x[i]-1)) + exp(-m * (T+1-x[i])) for i in 1:length(x)]
        isnothing(wpm) ? uwerr(cpp) : uwerr(cpp, wpm)                       
	    isnothing(wpm) ? uwerr.(R_dat) : [uwerr(R_dat[i], wpm) for i in 1:length(R_dat)]
        isnothing(wpm) ? uwerr.(vec) : [uwerr(vec[i], wpm) for i in 1:length(vec)]
        v = value.(vec)
        e = err.(vec)
        #v = value(cpp)
        #e = err(cpp)

        figure()
        errorbar(1:length(R_dat), value.(R_dat), err.(R_dat), fmt="x", color="black")
        fill_between(x, v-e, v+e, color="gray", alpha=0.75)
        ylabel(L"$R_\mathrm{PS}$")
        xlabel(L"$x_0$")
        #ylim(1e-8, 1e-7)
        #xlim(0,length(R_dat))
        yscale("log")
        tight_layout()
        savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.pdf"))

        #@gp x value.(R_dat) err.(R_dat) "w errorbars t ''"
        #@gp:- x v-e v+e "w filledcurves t ''"
        #@gp:- "set logscale y"
        #save(string("/home/asaez/cls_ens/codes/lattA.jl/plots/R_",ens.id,"_",PS,"_plat.gp"))
    end
    
    return f, syst, f_i, weight, pval    
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
        elseif ens.id in ["E250", "E300", "J500", "J501", "D450"]
            rwf = [read_ms1(path_rw[i], v=ens.vrw[i]) for i in 1:length(ens.vrw)]
            [Ysl[k] = Ysl[k][1:size(rwf[k],2), :, :] for k in 1:length(Ysl)]
        else
            rwf = read_ms1.(path_rw, v=ens.vrw)
        end

        Ysl_r, W = juobs.apply_rw(Ysl, rwf)
        tmp_r = Ysl_r[1]
        tmp_W = W[1]
        [tmp_r = cat(tmp_r, Ysl_r[k], dims=1) for k = 2:nr]
        [tmp_W = cat(tmp_W, W[k], dims=1) for k = 2:nr]
        W_obs = uwreal(tmp_W, id, replica, collect(1:length(tmp_W)), sum(replica))
        WY_aux = Matrix{uwreal}(undef, xmax, 2*dt0+1)
    end
    for i = 1:xmax
        k = 1
        for j = nt0-dt0:nt0+dt0
            if !rw
                Y_aux[i, k] = uwreal(tmp[:, i, j], id, replica)
            else
                WY_aux[i, k] = uwreal(tmp_r[:, i, j], id, replica, collect(1:length(tmp_W)), sum(replica))
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
                errorbar(1:length(dat), value.(dat), err.(dat), fmt="x", color="black")
                ylabel(L"$t^2\left<E\right>$")
                xlabel(L"$x_0$")
                ylim(v-5*e, v+20*e)
                tight_layout()

                savefig(string("/home/asaez/cls_ens/codes/lattA.jl/plots/t2E_",ens.id,"_plat.pdf"))
                
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
