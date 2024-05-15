using ADerrors, LsqFit, ForwardDiff, LinearAlgebra, SpecialFunctions, Optim

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

function fit_alg(f::Function, x::Union{Vector{Int64}, Vector{Float64}, Matrix{Float64}}, y::Vector{uwreal}, 
    n::Int64, guess::Union{Float64, Vector{Float64}, Nothing}=nothing; 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    isnothing(wpm) ? uwerr.(y) : [uwerr(y[i], wpm) for i in 1:length(y)]
    W = 1 ./ err.(y) .^ 2
    chisq = fit_defs(f,x,W)
    #lb = [-Inf for i in 1:n]
	#ub = [+Inf for i in 1:n]
    if guess == nothing
        p0 = [.5 for i in 1:n]
    else
        p0 = [guess; [1. for i in 1:n-length(guess)]]
        #lb[1] = .9 * guess  
        #ub[1] = 1.1 * guess
    end
    fit = curve_fit(f,x,value.(y),W,p0)#,lower=lb,upper=ub)
    chi2 = sum(fit.resid .^ 2)
    isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,coef(fit),y) : (up,chi_exp) = fit_error(chisq,coef(fit),y,wpm)
    isnothing(wpm) ? pval = pvalue(chisq,chi2,value.(up),y) : pval = pvalue(chisq,chi2,value.(up),y,wpm=wpm)
    return up, chi2, chi_exp, pval
end

function fit_alg(model::Function,x::Array{Float64},y::Array{uwreal},param::Int64,W::Matrix{Float64};
    guess::Union{Float64, Vector{Float64}, Nothing}=nothing)

    if guess == nothing
        p00 = [.5 for i in 1:param]
    else
        p00 = [guess; [1. for i in 1:param-length(guess)]]
        #lb[1] = .9 * guess  
        #ub[1] = 1.1 * guess
    end

    chisq_corr(par,dat) = juobs.gen_chisq_correlated(model, x, W, par, dat)
    fit = curve_fit(model,x,value.(y),W,p00)
    up, chi_exp = fit_error(chisq_corr, coef(fit), y, W=W)
    uwerr.(up)
    chi2 = sum(fit.resid .^ 2)
    pval = pvalue(chisq_corr, sum(fit.resid .^ 2), value.(up), y, W)
    doff = dof(fit)

    return up, chi2, chi_exp, pval, doff
end

#=
function fit_alg(model::Function,x::Array{Float64},y::Array{uwreal},param::Int64,W::Matrix{Float64})
    p00 = [0.5 for i in 1:param]

    chisq_corr(par,dat) = juobs.gen_chisq_correlated(model, x, W, par, dat)
    f(p) = LsqFit.cholesky(W).U * (model(x, p) - value.(y)) 
    fit = LsqFit.lmfit(f,p00,W) 
    up, chi_exp = fit_error(chisq_corr, coef(fit), y, W=W)
    uwerr.(up)
    chi2 = sum(fit.resid .^ 2)
    pval = pvalue(chisq_corr, sum(fit.resid .^ 2), value.(up), y, W)
    doff = dof(fit)

    return up, chi2, chi_exp, pval, doff
end
=#

function fit_alg(f::Vector{Function}, x::Vector{Matrix{Float64}}, y::Vector{Vector{uwreal}}, 
    n::Int64, guess::Union{Float64, Vector{Float64}, Nothing}=nothing; 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    y_n = y[1]
    x_n = x[1]
    for i in 2:length(x)
        y_n = vcat(y_n, y[i])
        x_n = vcat(x_n, x[i])
    end
    idx = [collect(1:size(x[i],1)) for i in 1:length(x)]
    [idx[i] = idx[i] .+ idx[i-1][end] for i in 2:length(idx)]

    isnothing(wpm) ? [uwerr.(y[i]) for i in 1:length(y)] : [[uwerr(y[i][j], wpm) for j in 1:length(y[i])] for i in 1:length(y)]
    W = [1 ./ err.(y[i]) .^ 2 for i in 1:length(y)]

    chisq = (par, y_n) -> sum([sum((y_n[idx[i]] .- f[i](x[i], par)) .^ 2 .* W[i]) for i in 1:length(x)])
    min_fun(t) = chisq(t, value.(y_n))
    if guess == nothing
        p0 = [.5 for i in 1:n]
    else
        p0 = [guess; [1. for i in 1:n-length(guess)]]
    end
    sol = optimize(min_fun, p0, Optim.Options(g_tol=1e-8, iterations=1000000)) 
    chi2 = min_fun(sol.minimizer)
    isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n) : (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n,wpm)
    isnothing(wpm) ? pval = pvalue(chisq,chi2,value.(up),y_n) : pval = pvalue(chisq,chi2,value.(up),y_n,wpm=wpm)
    return up, chi2, chi_exp, pval
end

##CHECK THIS ONE, NOT WORKING GOOD
function fit_alg(f::Vector{Function}, x::Vector{Matrix{Float64}}, y::Vector{Vector{uwreal}}, W::Vector{Matrix{Float64}}, 
    n::Int64, guess::Union{Float64, Vector{Float64}, Nothing}=nothing; 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    y_n = y[1]
    x_n = x[1]
    for i in 2:length(x)
        y_n = vcat(y_n, y[i])
        x_n = vcat(x_n, x[i])
    end
    idx = [collect(1:size(x[i],1)) for i in 1:length(x)]
    [idx[i] = idx[i] .+ idx[i-1][end] for i in 2:length(idx)]

    isnothing(wpm) ? [uwerr.(y[i]) for i in 1:length(y)] : [[uwerr(y[i][j], wpm) for j in 1:length(y[i])] for i in 1:length(y)]

    chisq = (par, y_n) -> sum([juobs.gen_chisq_correlated(f[i], x[i], W[i], par, y_n[i]) for i in 1:length(x)])

    min_fun(t) = chisq(t, value.(y_n))
    if guess == nothing
        p0 = [.5 for i in 1:n]
    else
        p0 = [guess; [1. for i in 1:n-length(guess)]]
    end
    sol = optimize(min_fun, p0, Optim.Options(g_tol=1e-8, iterations=1000000)) 
    chi2 = min_fun(sol.minimizer)

    W_aux = Matrix{Float64}(undef, length(y_n), length(y_n))
    for i in 1:length(y_n)
        for j in 1:length(y_n)
            if i <= length(y[1])
                if j <= length(y[1])
                    W_aux[i,j] = W[1][i,j]
                else
                    W_aux[i,j] = .0
                end
            else
                if j <= length(y[1])
                    W_aux[i,j] = .0
                else
                    W_aux[i,j] = W[2][i-length(y[1]),j-length(y[1])]
                end
            end
        end
    end

    isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n,W=W_aux) : (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n,wpm,W=W_aux)
    isnothing(wpm) ? pval = pvalue(chisq,chi2,value.(up),y_n,W_aux) : pval = pvalue(chisq,chi2,value.(up),y_n,wpm=wpm,W_aux)
    return up, chi2, chi_exp, pval
end

function fit_alg(f::Vector{Function}, x::Vector{Matrix{Float64}}, y::Vector{Vector{uwreal}}, 
    n::Int64, lb::Vector{Float64}, ub::Vector{Float64}, guess::Union{Float64, Vector{Float64}, Nothing}=nothing; 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    y_n = y[1]
    x_n = x[1]
    for i in 2:length(x)
        y_n = vcat(y_n, y[i])
        x_n = vcat(x_n, x[i])
    end
    idx = [collect(1:size(x[i],1)) for i in 1:length(x)]
    [idx[i] = idx[i] .+ idx[i-1][end] for i in 2:length(idx)]

    isnothing(wpm) ? [uwerr.(y[i]) for i in 1:length(y)] : [[uwerr(y[i][j], wpm) for j in 1:length(y[i])] for i in 1:length(y)]
    W = [1 ./ err.(y[i]) .^ 2 for i in 1:length(y)]

    chisq = (par, y_n) -> sum([sum((y_n[idx[i]] .- f[i](x[i], par)) .^ 2 .* W[i]) for i in 1:length(x)])
    min_fun(t) = chisq(t, value.(y_n))
    if guess == nothing
        p0 = [.5 for i in 1:n]
    else
        p0 = [guess; [1. for i in 1:n-length(guess)]]
    end
    sol = optimize(min_fun, lb, ub, p0, Fminbox(NelderMead()), Optim.Options(g_tol=1e-8, iterations=1000000))
    chi2 = min_fun(sol.minimizer)
    isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n) : (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n,wpm)
    isnothing(wpm) ? pval = pvalue(chisq,chi2,value.(up),y_n) : pval = pvalue(chisq,chi2,value.(up),y_n,wpm=wpm)
    return up, chi2, chi_exp, pval
end

function fit_alg_LBFGS(f::Vector{Function}, x::Vector{Matrix{Float64}}, y::Vector{Vector{uwreal}}, 
    n::Int64, guess::Union{Float64, Vector{Float64}, Nothing}=nothing; 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    
    y_n = y[1]
    x_n = x[1]
    for i in 2:length(x)
        y_n = vcat(y_n, y[i])
        x_n = vcat(x_n, x[i])
    end
    idx = [collect(1:size(x[i],1)) for i in 1:length(x)]
    [idx[i] = idx[i] .+ idx[i-1][end] for i in 2:length(idx)]

    isnothing(wpm) ? [uwerr.(y[i]) for i in 1:length(y)] : [[uwerr(y[i][j], wpm) for j in 1:length(y[i])] for i in 1:length(y)]
    W = [1 ./ err.(y[i]) .^ 2 for i in 1:length(y)]

    chisq = (par, y_n) -> sum([sum((y_n[idx[i]] .- f[i](x[i], par)) .^ 2 .* W[i]) for i in 1:length(x)])
    min_fun(t) = chisq(t, value.(y_n))
    if guess == nothing
        p0 = [.5 for i in 1:n]
    else
        p0 = [guess; [1. for i in 1:n-length(guess)]]
    end
    sol = optimize(min_fun, p0, LBFGS()) 
    chi2 = min_fun(sol.minimizer)
    isnothing(wpm) ? (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n) : (up,chi_exp) = fit_error(chisq,sol.minimizer,y_n,wpm)
    isnothing(wpm) ? pval = pvalue(chisq,chi2,value.(up),y_n) : pval = pvalue(chisq,chi2,value.(up),y_n,wpm=wpm)
    return up, chi2, chi_exp, pval
end

function model_av(fun::Function, y::Vector{uwreal}, guess::Float64; 
    tm::Vector{Int64}, tM::Vector{Int64}, k::Int64, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing) 
    
    pval = Array{Float64,1}()
    p_1 = Array{uwreal,1}()
    TIC = Array{Float64,1}()

    isnothing(wpm) ? uwerr.(y) : [uwerr(y[i], wpm) for i in 1:length(y)]
    for i in tm
        for j in tM
            if i < j
                x = collect(i:j)
                y_aux = y[i:j]
                try 
                    up, chi2, chi_exp, pval_i = fit_alg(fun,x,y_aux,k,guess,wpm=wpm)
                    push!(pval, pval_i)
                    push!(TIC, chi2 - 2*chi_exp)
                    push!(p_1, up[1])
                catch e 
                end
            end
        end
    end

    TIC = TIC .- minimum(TIC)
    weight = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
    p_av = sum(p_1 .* weight)
    syst2 = sum(p_1 .^ 2 .* weight) - p_av ^ 2
    return p_av, syst2, p_1, weight, pval
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
                if i < j
                    x = collect(i:j)
                    y_aux = y[i:j]
                    try 
                        up, chi2, chi_exp, pval_i = fit_alg(f,x,y_aux,k[ind],guess,wpm=wpm)
                        push!(pval, pval_i)
                        push!(TIC, chi2 - 2*chi_exp)
                        push!(p_1, up[1])
                    catch e 
                    end
                end
            end
        end
    end

    TIC = TIC .- minimum(TIC)
    weight = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
    p_av = sum(p_1 .* weight)
    syst2 = sum(p_1 .^ 2 .* weight) - p_av ^ 2
    return p_av, syst2, p_1, weight, pval
end

function pvalue(chisq::Function,
    chi2::Float64,
    xp::Vector{Float64}, 
    data::Vector{uwreal};
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} =  Dict{Int64,Vector{Float64}}(),
    nmc::Int64 = 5000)

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

    if (m-n > 0)
        W = zeros(Float64, m)
            for i in 1:m
                if (data[i].err == 0.0)
                    #isnothing(wpm) ? wuerr(data[i]) : uwerr(data[i], wpm)
                    uwerr(data[i], wpm)
                    if (data[i].err == 0.0)
                        error("Zero error in fit data")
                    end
                end
                global W[i] = 1.0 / data[i].err^2
            end

        m = length(data)
        n = size(hess, 1) - m

        hm = view(hess, 1:n, n+1:n+m)
        sm = Array{Float64, 2}(undef, n, m)
        for i in 1:n, j in 1:m
            sm[i,j] = hm[i,j] / sqrt.(W[j])
        end
        maux = sm * sm'
        hi   = LinearAlgebra.pinv(maux)
        Px   = -hm' * hi * hm

        for i in 1:m
            Px[i,i] = W[i] + Px[i,i]
        end

        C = cov(data) 
        
        nu = sqrt(C) * Px * sqrt(C)
        
        N = length(nu[1,:])
        z = randn(N, nmc)

        eig = abs.(eigvals(nu))
        eps = 1e-14 * maximum(eig)
        eig = eig .* (eig .> eps)

        aux = eig' * (z .^ 2)
        global Q = 1.0 - juobs.mean(aux .< chi2)

        x = chi2 .- eig[2:end]' * (z[2:end,:].^2)
        x = x / eig[1]
        #dQ = juobs.mean((x .> 0) .* exp.(-x * 0.5) * 0.5 ./ sqrt.(abs.(x)))
        #dQ = err(cse)/value(cse) * dQ
    end
    return Q
end

function pvalue(chisq::Function,
    chi2::Float64,
    xp::Vector{Float64}, 
    data::Vector{uwreal},
    W::Array{Float64,2};
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing} =  Dict{Int64,Vector{Float64}}(),
    nmc::Int64 = 5000)

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

    if (m-n > 0)
        m = length(data)
        n = size(hess, 1) - m
        
        Lm = LinearAlgebra.cholesky(LinearAlgebra.Symmetric(W))
        Li = LinearAlgebra.inv(Lm.L)
        
        hm = view(hess, 1:n, n+1:n+m)
        sm = hm * Li'
        
        maux = sm * sm'
        hi   = LinearAlgebra.pinv(maux)
        Px   = W - hm' * hi * hm
        
        C = cov(data) 
        
        nu = sqrt(C) * Px * sqrt(C)
        
        N = length(nu[1,:])
        z = randn(N, nmc)

        eig = abs.(eigvals(nu))
        eps = 1e-14 * maximum(eig)
        eig = eig .* (eig .> eps)

        aux = eig' * (z .^ 2)
        global Q = 1.0 - juobs.mean(aux .< chi2)

        x = chi2 .- eig[2:end]' * (z[2:end,:].^2)
        x = x / eig[1]
        #dQ = juobs.mean((x .> 0) .* exp.(-x * 0.5) * 0.5 ./ sqrt.(abs.(x)))
        #dQ = err(cse)/value(cse) * dQ
    end
    return Q
end


function fve(mpi::uwreal, mk::uwreal, fpi::uwreal, fk::uwreal, ens::EnsInfo)
    mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    meta = sqrt(4/3 * mk ^ 2 - 1/3 * mpi ^ 2)

    jipi = value(mpi) ^ 2 / (4*pi*fpi) ^ 2
    jik = value(mk) ^ 2 / (4*pi*fpi) ^ 2
    jieta = value(meta) ^ 2 / (4*pi*fpi) ^ 2
    mpiL = value(mpi) * ens.L
    mkL = value(mk) * ens.L
    metaL = value(meta) * ens.L
    lampi = mpiL * sqrt.(nn)
    lamk = mkL * sqrt.(nn)
    lameta = metaL * sqrt.(nn)
    g1pi = sum(4 .* mm ./ lampi .* besselk.(1, lampi))
    g1k = sum(4 .* mm ./ lamk .* besselk.(1, lamk))
    g1eta = sum(4 .* mm ./ lameta .* besselk.(1, lameta))
    fve_mpi = 1/4 * jipi * g1pi - 1/12 * jieta * g1eta
    fve_mk = 1/6 * jieta * g1eta
    fve_fpi = - jipi * g1pi - 1/2 * jik * g1k
    fve_fk = -3/8 * jipi *g1pi - 3/4 * jik * g1k - 3/8 * jieta * g1eta

    mpi_infty = mpi / (1+fve_mpi)
    mk_infty = mk / (1+fve_mk)
    fpi_infty = fpi / (1+fve_fpi)
    fk_infty = fk / (1+fve_fk)

    return mpi_infty, mk_infty, fpi_infty, fk_infty
end

function fve_inv(mpi::uwreal, mk::uwreal, fpi::uwreal, fk::uwreal, ens::EnsInfo)
    mm = [6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24]
	nn = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    meta = sqrt(4/3 * mk ^ 2 - 1/3 * mpi ^ 2)

    jipi = value(mpi) ^ 2 / (4*pi*fpi) ^ 2
    jik = value(mk) ^ 2 / (4*pi*fpi) ^ 2
    jieta = value(meta) ^ 2 / (4*pi*fpi) ^ 2
    mpiL = value(mpi) * ens.L
    mkL = value(mk) * ens.L
    metaL = value(meta) * ens.L
    lampi = mpiL * sqrt.(nn)
    lamk = mkL * sqrt.(nn)
    lameta = metaL * sqrt.(nn)
    g1pi = sum(4 .* mm ./ lampi .* besselk.(1, lampi))
    g1k = sum(4 .* mm ./ lamk .* besselk.(1, lamk))
    g1eta = sum(4 .* mm ./ lameta .* besselk.(1, lameta))
    fve_mpi = 1/4 * jipi * g1pi - 1/12 * jieta * g1eta
    fve_mk = 1/6 * jieta * g1eta
    fve_fpi = -1 * jipi * g1pi - 1/2 * jik * g1k
    fve_fk = -3/8 * jipi *g1pi - 3/4 * jik * g1k - 3/8 * jieta * g1eta

    mpi_L = mpi * (1+fve_mpi)
    mk_L = mk * (1+fve_mk)
    fpi_L = fpi * (1+fve_fpi)
    fk_L = fk * (1+fve_fk)

    return mpi_L, mk_L, fpi_L, fk_L
end

function corr_sym_E250(corr1::juobs.Corr, corr2::juobs.Corr, parity::Int64=1)
    aux = [corr2.obs[97:end]; corr2.obs[1:96]]
    corr2_sym = juobs.Corr(aux, corr2.kappa, corr2.mu, corr2.gamma, corr1.y0, corr2.theta1, corr2.theta2)

    corr = [corr1.obs[1:3]; (corr1.obs[4:98] .+ parity * corr2_sym.obs[192:-1:98]) / 2]

    return juobs.Corr(corr, corr1.kappa, corr1.mu, corr1.gamma, corr1.y0, corr1.theta1, corr1.theta2)
end

function corr_sym_D450(corr1::juobs.Corr, corr2::juobs.Corr, parity::Int64=1)
    aux = [corr2.obs[65:end]; corr2.obs[1:64]]
    corr2_sym = juobs.Corr(aux, corr2.kappa, corr2.mu, corr2.gamma, corr1.y0, corr2.theta1, corr2.theta2)

    corr = [corr1.obs[1:3]; (corr1.obs[4:66] .+ parity * corr2_sym.obs[128:-1:66]) / 2]

    return juobs.Corr(corr, corr1.kappa, corr1.mu, corr1.gamma, corr1.y0, corr1.theta1, corr1.theta2)
end
