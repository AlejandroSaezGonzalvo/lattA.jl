using ADerrors, LsqFit

function fit_defs(f::Function,x,W) ## uncorrelated fit
	chisq(p,d) = sum((d .- f(x,p)).^2 .* W)
	return chisq
end

function model_av(f::Function, y::Vector{uwreal}, guess::Float64; tm::Vector{Int64}, tM::Vector{Int64}, k::Int64, wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing) 
    pval = Array{Float64,1}()
    p_1 = Array{uwreal,1}()
    TIC = Array{Float64,1}()

    isnothing(wpm) ? uwerr.(y) : [uwerr(y[i], wpm) for i in 1:length(y)]
    for i in tm
        for j in tM
            x = collect(i:j)
            dy = err.(y)
            W = 1 ./ dy .^ 2
            p0 = [guess; [0.5 for i in 1:k-1]]
            chisq = fit_defs(f,x,W)
            fit = curve_fit(f,x,value.(y),W,p0)
            chi2 = sum(fit.resid .^ 2)
            isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,coef(fit),y) : (up,chi_exp) = fit_error(chisq,coef(fit),y,wpm)
            isnothing(wpm) ? push!(pval, pvalue(chisq,chi2,value.(up),y)) : push!(pval, pvalue(chisq,chi2,value.(up),y,wpm=wpm))
            push!(TIC, chi2 - 2*chi_exp)
            push!(p_1, up[1])
        end
    end

    TIC = TIC .- minimum(TIC)
    weight = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
    p_av = sum(p_1 .* weight)
    syst2 = sum(p_1 .^ 2 .* weight) .- p_av .^ 2
    return p_av, syst2, p_1, weight, pval
end
