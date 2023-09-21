using juobs

function get_m(corr::Corr, id::String; tm::Union{Vector{Int64}, Nothing}=nothing, tM::Union{Vector{Int64}, Nothing}=nothing, pl::Bool=false, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    corr_d = corr.obs
    m_dat = 0.5 .* log.((corr_d[1:end-1] ./ corr_d[2:end]) .^ 2)
    y0 = corr.y0
    T = length(corr_d) - y0
    
    isnothing(tm) ? tm = [y0+10, y0+15, y0+20, y0+25] : tm=tm
    isnothing(tM) ? tM = [T-10, T-15, T-20, T-25] : tM=tM
    @.fit_fun(x,p) = p[1] + p[2] * exp(-p[3] * (x-y0)) + p[4] * exp(-p[5] * (T-x))
    k = 5
    
    m, syst2, m_i, weight, pval = model_av(fit_fun, m_dat, ens_obs[id][2], tm=tm, tM=tM, k=k, wpm=wpm)
    if pl == true
        ##plot
        savefig(string("/home/asaez/cls_ens/codes/analysis_cls/plots/m_",id,"_",round(m,digits=3),".pdf"))
        close("all")
    end

    return m, syst2
end