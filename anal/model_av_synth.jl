using juobs, ADerrors, LsqFit

function fit_defs(f::Function,x,W) 
	chisq(p,d) = sum((d .- f(x,p)).^2 .* W)
	return chisq
end

function fit_alg(model::Function, xdata::Array{<:Real}, ydata::Array{uwreal}, param::Int64=3; wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing)
    isnothing(wpm) ? uwerr.(ydata) : uwerr.(ydata, wpm)
    
    yval = value.(ydata)
    yer = err.(ydata)
    W = 1 ./ yer .^ 2
    
    # Generate chi2 + solver
    chisq = fit_defs(model, xdata, W)
    if param == length(ptrue)
        pp = ptrue
    else
        pp = [ptrue; fill(0.5, param-length(ptrue))]
    end
    fit = curve_fit(model, xdata, yval, 1.0 ./ yer.^2, pp)
    upar, chi_exp = isnothing(wpm) ? fit_error(chisq, coef(fit), ydata) : fit_error(chisq, coef(fit), ydata, wpm)
    chi2 = sum(fit.resid .^ 2)
    pval = pvalue(chisq, sum(fit.resid .^ 2), value.(upar), ydata)
    doff = dof(fit)
    
    return upar, chi_exp, chi2, pval, doff 
end

function fpi_SU2(x,p) 
    return [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,1] for i in 1:length(x[:,2])]
end

function fpi_SU2_a2phi2(x,p) 
    return [p[1] * x[i,2] + p[2] / (4 * pi) * (1 - 2 * x[i,2] / p[2] ^ 2 * log(x[i,2] / p[2] ^ 2)) + p[3] * x[i,1] + p[4] * x[i,1] * x[i,2] for i in 1:length(x[:,2])]
end

################################################################

phi2beta1 = collect(0.1:0.05:0.67)
phi2beta2 = collect(0.07:0.05:0.67)
phi2beta3 = collect(0.12:0.05:0.67)
phi2beta4 = collect(0.11:0.05:0.67)
phi2beta5 = collect(0.23:0.05:0.67)
t0 = collect(2.8:2.5:13)
x = [[[1 / t0[1] for i in 1:length(phi2beta1)] phi2beta1]; 
     [[1 / t0[2] for i in 1:length(phi2beta2)] phi2beta2];
     [[1 / t0[3] for i in 1:length(phi2beta3)] phi2beta3];
     [[1 / t0[4] for i in 1:length(phi2beta4)] phi2beta4];
     [[1 / t0[5] for i in 1:length(phi2beta5)] phi2beta5]]

ptrue = [-0.04,3.2,-0.2]
xph = [0.0 0.1]

n = length(x[:,1])
eta = [uwreal(randn(10000), i+1) for i in 1:n]
dat = fpi_SU2(x,ptrue) .+ eta
uwerr.(dat)

cut = [collect(1:n), collect(1:2:n), collect(1:3:n), collect(1:4:n), [1,5,6,7,10,20,23,24,30,35,40,41,42,43,46,47,50,51,52,53,54,55,56,57,58], [1,3,5,8,15,20,22,25], collect(50:58), collect(40:58)]

y = [dat[cut[i]] for i in 1:length(cut)]
xcut = [x[cut[i],:] for i in 1:length(cut)]

models = [fpi_SU2, fpi_SU2_a2phi2]
params = [3,4]

TIC = Array{Float64,1}()
dof_vec = Array{Float64,1}()
yph = Array{uwreal,1}()

for j in 1:length(cut)
    for i in 1:length(models)
        upar, chi_exp, chi2, pval, doff = fit_alg(models[i], xcut[j], y[j], params[i])
        push!(TIC, chi2 - 2 * chi_exp)
        push!(dof_vec, doff)
        push!(yph, models[i](xph, upar)[1])
        println(i, " ", j)
    end
end
W = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
yph_av = sum(yph .* W)
syst2 = sum(yph .^ 2 .* W) - yph_av ^ 2

#########################################################################

phi2beta1 = collect(0.2:0.05:0.67)
phi2beta2 = collect(0.17:0.05:0.67)
phi2beta3 = collect(0.21:0.05:0.67)
phi2beta4 = collect(0.22:0.05:0.67)
t0 = collect(2.8:2.5:13)
x = [[[1 / t0[1] for i in 1:length(phi2beta1)] phi2beta1]; 
     [[1 / t0[2] for i in 1:length(phi2beta2)] phi2beta2];
     [[1 / t0[3] for i in 1:length(phi2beta3)] phi2beta3];
     [[1 / t0[4] for i in 1:length(phi2beta4)] phi2beta4]]

ptrue = [-0.04,3.2,-0.2]
xph = [0.0 0.1]

n = length(x[:,1])
eta = [uwreal(randn(10000), i+1) for i in 1:n]
dat = fpi_SU2(x,ptrue) .+ eta
uwerr.(dat)

cut = [collect(1:n), collect(1:2:n), collect(1:3:n), collect(1:4:n), [1,5,6,7,10,20,23,24,30,35,40,41], [1,3,5,8,15,20,22,25]]

y = [dat[cut[i]] for i in 1:length(cut)]
xcut = [x[cut[i],:] for i in 1:length(cut)]

models = [fpi_SU2, fpi_SU2_a2phi2]
params = [3,4]

TIC = Array{Float64,1}()
dof_vec = Array{Float64,1}()
yph = Array{uwreal,1}()

for j in 1:length(cut)
    for i in 1:length(models)
        upar, chi_exp, chi2, pval, doff = fit_alg(models[i], xcut[j], y[j], params[i])
        push!(TIC, chi2 - 2 * chi_exp)
        push!(dof_vec, doff)
        push!(yph, models[i](xph, upar)[1])
        println(i, " ", j)
    end
end
W = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
yph2_av = sum(yph .* W)
syst22 = sum(yph .^ 2 .* W) - yph_av ^ 2

################################################################

phi2beta1 = collect(0.1:0.05:0.67)
phi2beta2 = collect(0.07:0.05:0.67)
phi2beta3 = collect(0.12:0.05:0.67)
phi2beta4 = collect(0.11:0.05:0.67)
phi2beta5 = collect(0.23:0.05:0.67)
t0 = collect(2.8:2.5:13)
x = [[[1 / t0[1] for i in 1:length(phi2beta1)] phi2beta1]; 
     [[1 / t0[2] for i in 1:length(phi2beta2)] phi2beta2];
     [[1 / t0[3] for i in 1:length(phi2beta3)] phi2beta3];
     [[1 / t0[4] for i in 1:length(phi2beta4)] phi2beta4];
     [[1 / t0[5] for i in 1:length(phi2beta5)] phi2beta5]]

ptrue = [-0.04,3.2,-0.2]
xph = [0.0 0.1]

n = length(x[:,1])
eta = [uwreal(randn(10000), i+1) for i in 1:n]
dat = fpi_SU2(x,ptrue) .+ eta
uwerr.(dat)

cut = [collect(1:3:n), collect(1:4:n), [1,5,6,7,10,20,23,24,30,35,40,41,42,43,46,47,50,51,52,53,54,55,56,57,58], [1,3,5,8,15,20,22,25], collect(50:58), collect(40:58)]

y = [dat[cut[i]] for i in 1:length(cut)]
xcut = [x[cut[i],:] for i in 1:length(cut)]

models = [fpi_SU2, fpi_SU2_a2phi2]
params = [3,4]

TIC = Array{Float64,1}()
dof_vec = Array{Float64,1}()
yph = Array{uwreal,1}()

for j in 1:length(cut)
    for i in 1:length(models)
        upar, chi_exp, chi2, pval, doff = fit_alg(models[i], xcut[j], y[j], params[i])
        push!(TIC, chi2 - 2 * chi_exp)
        push!(dof_vec, doff)
        push!(yph, models[i](xph, upar)[1])
        println(i, " ", j)
    end
end
W = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC))
yph3_av = sum(yph .* W)
syst32 = sum(yph .^ 2 .* W) - yph_av ^ 2
