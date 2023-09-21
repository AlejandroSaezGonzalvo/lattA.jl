using juobs, LaTeXStrings, PyPlot
include("const.jl");

function get_m(corr::juobs.Corr, id::String; tm::Union{Vector{Int64}, Nothing}=nothing, tM::Union{Vector{Int64}, Nothing}=nothing, pl::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    corr_d = corr.obs
    m_dat = 0.5 .* log.((corr_d[2:end-2] ./ corr_d[3:end-1]) .^ 2)
    y0 = corr.y0
    T = length(corr_d) - y0
    
    isnothing(tm) ? tm = [y0+10,y0+12,y0+14] : tm=tm
    isnothing(tM) ? tM = [T-10,T-12,T-14] : tM=tM
    @.fit_fun(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-x))
    k = 5
    
    m, syst, m_i, weight, pval = model_av(fit_fun, m_dat, ens_obs[id][2], tm=tm, tM=tM, k=k, wpm=wpm)
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
        ylim(v-50*e, v+100*e)

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

        savefig(string("/home/asaez/cls_ens/codes/analysis_cls/plots/m_",id,"_",round(value(m),digits=3),".pdf"))
        close("all")
    end

    return m, syst, m_i, weight, pval
end
