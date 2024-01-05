using ADerrors, LsqFit, ForwardDiff, LinearAlgebra, SpecialFunctions

function get_model(x, p, n)
    s = 0.0
    for k = 1:n
        s = s .+ p[k] .* x.^(k-1)
    end
    return s
end

function fit_defs(f::Function,x,W) ## uncorrelated fit
	chisq(p,d) = sum((d .- f(x,p)).^2 .* W)
	return chisq
end

function fit_alg(f::Function, x::Union{Vector{Int64}, Vector{Float64}}, y::Vector{uwreal}, 
    n::Int64, guess::Union{Float64, Nothing}=nothing; 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    isnothing(wpm) ? uwerr.(y) : [uwerr(y[i], wpm) for i in 1:length(y)]
    W = 1 ./ err.(y) .^ 2
    chisq = fit_defs(f,x,W)
    #lb = [-Inf for i in 1:n]
	#ub = [+Inf for i in 1:n]
    if guess == nothing
        p0 = [.5 for i in 1:n]
    else
        p0 = [guess; [1. for i in 1:n-1]]
        #lb[1] = .9 * guess  
        #ub[1] = 1.1 * guess
    end
    fit = curve_fit(f,x,value.(y),W,p0)#,lower=lb,upper=ub)
    chi2 = sum(fit.resid .^ 2)
    isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,coef(fit),y) : (up,chi_exp) = fit_error(chisq,coef(fit),y,wpm)
    isnothing(wpm) ? pval = pvalue(chisq,chi2,value.(up),y) : pval = pvalue(chisq,chi2,value.(up),y,wpm=wpm)
    return up, chi2, chi_exp, pval
end

function model_av(fun::Vector{Function}, y::Vector{uwreal}, guess::Float64; 
    tm::Vector{Vector{Int64}}, tM::Vector{Vector{Int64}}, k::Vector{Int64}, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing) 
    
    pval = Array{Float64,1}()
    p_1 = Array{uwreal,1}()
    TIC = Array{Float64,1}()

    isnothing(wpm) ? uwerr.(y) : [uwerr(y[i], wpm) for i in 1:length(y)]
    for ind in 1:length(fun)
        f = fun[ind]
        for i in tm[ind]
            for j in tM[ind]
                x = collect(i:j)
                y_aux = y[i:j]
                try 
                    up, chi2, chi_exp, pval_i = fit_alg(f,x,y_aux,k[ind],guess,wpm=wpm)
                catch e 
                end
                push!(pval, pval_i)
                push!(TIC, chi2 - 2*chi_exp)
                push!(p_1, up[1])
            end
        end
    end

    TIC = TIC .- minimum(TIC)
    weight = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
    p_av = sum(p_1 .* weight)
    syst = sqrt(sum(p_1 .^ 2 .* weight) - p_av ^ 2)
    return p_av, syst, p_1, weight, pval
end

function pvalue(chisq::Function, chi2::Float64, xp::Vector{Float64}, data::Vector{uwreal};
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=Dict{String,Vector{Float64}}(),
    W::Union{Vector{Float64},Array{Float64,2}} = Vector{Float64}(), nmc::Int64=5000)

    n = length(xp)   # Number of fit parameters
    m = length(data) # Number of data

    xav = zeros(Float64, n+m)
    for i in 1:n
        xav[i] = xp[i]
    end
    for i in n+1:n+m
        xav[i] = data[i-n].mean
    end
    ccsq(x::Vector) = chisq(view(x, 1:n), view(x, n+1:n+m)) 
    if (n+m < 4)
        cfg = ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{1}());
    else
        cfg = ForwardDiff.HessianConfig(ccsq, xav, ADerrors.Chunk{4}());
    end
        
    hess = Array{Float64}(undef, n+m, n+m)
    ForwardDiff.hessian!(hess, ccsq, xav, cfg)
        
    cse = 0.0
    Q = dQ = 0.0
    if (m-n > 0)
        if (length(W) == 0)
            Ww = zeros(Float64, m)
            for i in 1:m
                if (data[i].err == 0.0)
                    uwerr(data[i], wpm)
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                Ww[i] = 1.0 / data[i].err^2
            end
        else
            Ww = W
        end
        #cse = chiexp(hess, data, Ww, wpm)

        m = length(data)
        n = size(hess, 1) - m

        hm = view(hess, 1:n, n+1:n+m)
        sm = Array{Float64, 2}(undef, n, m)

        for i in 1:n, j in 1:m
            sm[i,j] = hm[i,j] / sqrt.(Ww[j])
        end
        maux = sm * sm'
        hi   = LinearAlgebra.pinv(maux)
        Px   = -hm' * hi * hm

        for i in 1:m
            Px[i,i] = Ww[i] + Px[i,i]
        end
        C = cov(data, wpm) 
        
        nu = sqrt(C) * Px * sqrt(C)
        
        N = length(nu[1,:])
        z = randn(N, nmc)

        eig = abs.(eigvals(nu))
        eps = 1e-14 * maximum(eig)
        eig = eig .* (eig .> eps)

        aux = eig' * (z .^ 2)
        Q = 1.0 - juobs.mean(aux .< chi2)

        x = chi2 .- eig[2:end]' * (z[2:end,:].^2)
        x = x / eig[1]
        #dQ = juobs.mean((x .> 0) .* exp.(-x * 0.5) * 0.5 ./ sqrt.(abs.(x)))
        #dQ = err(cse)/value(cse) * dQ
    end

    return Q #uwreal([Q,dQ],"")
end

function fve(mpi::uwreal, mk::uwreal, fpi::uwreal, fk::uwreal, ens::EnsInfo)

    mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    jipi = value(mpi) / (4*pi*fpi) ^ 2
    jik = value(mk) / (4*pi*fk) ^ 2
    mpiL = value(mpi) * ens.L
    mkL = value(mk) * ens.L
    lampi = mpiL * sqrt.(nn)
    lamk = mkL * sqrt.(nn)
    g1pi = sum(4 .* mm ./ lampi .* besselk.(1, lampi))
    g1k = sum(4 .* mm ./ lamk .* besselk.(1, lamk))
    fve_mpi = 0.5 * jipi * g1pi
    fve_fpi = -2 * jipi * g1pi - jik * g1k
    fve_fk = -3/4 * jipi *g1pi - 3/2 * jik * g1k

    mpi = mpi / (1+fve_mpi)
    fpi = fpi / (1+fve_fpi)
    fk = fk / (1+fve_fk)

    return mpi, fpi, fk
end
