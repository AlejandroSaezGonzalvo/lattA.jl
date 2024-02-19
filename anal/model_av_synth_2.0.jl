using juobs, ADerrors, LsqFit, PyPlot

hc = 197.3269804 ## 1 MeV * fm = 1/hc
Lambda = 0.4841473 / hc * 1000
Nf = 3
b0 = (33-2*Nf)/(12*pi)
LambdaMS = 0.341 / hc * 1000
G1 = 0.24731
G2 = 1.63762 

fpi_exp = uwreal([130.56,0.02],"fpi_exp") + uwreal([0.0,0.13],"fpi QEDD") + uwreal([0.0,0.02],"fpi VUD") ## MeV ; FLAG
fk_exp = uwreal([157.2,0.2],"fk_expp") + uwreal([0.0,0.2],"fk QEDD") + uwreal([0.0,0.4],"fk VUD") ## MeV
fpik_exp = (2/3) * ( 0.5 * fpi_exp + fk_exp )
fpik_exp = fpik_exp / hc
uwerr(fpik_exp)

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
    elseif param > length(ptrue) 
        pp = [ptrue; fill(0.5, param-length(ptrue))]
    else 
        pp = ptrue[1:param]
    end
    fit = curve_fit(model, xdata, yval, 1.0 ./ yer.^2, pp)
    upar, chi_exp = isnothing(wpm) ? fit_error(chisq, coef(fit), ydata) : fit_error(chisq, coef(fit), ydata, wpm)
    chi2 = sum(fit.resid .^ 2)
    pval = pvalue(chisq, sum(fit.resid .^ 2), value.(upar), ydata)
    doff = dof(fit)
    
    return upar, chi_exp, chi2, pval, doff 
end

function SU3_continuum(x,p)
    f = [(p[1] / (4 * pi)) * (1 - 7/6 * (x[i,2]/p[1]^2*log(x[i,2]/p[1]^2)) - 4/3 * ((x[i,3]-1/2*x[i,2])/p[1]^2*log((x[i,3]-1/2*x[i,2])/p[1]^2)) - 1/2 * ((4/3*x[i,3]-x[i,2])/p[1]^2*log((4/3*x[i,3]-x[i,2])/p[1]^2)) + p[2] * x[i,3]) for i in 1:length(x[:,1])]
    return f
end

function Tay_continuum(x,p)
    f = [p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2 for i in 1:length(x[:,1])]
    return f
end

function Tay4_continuum(x,p)
    f = [p[1] + p[2] * (x[i,2] - x[i,4]) ^ 2 + p[3] * (x[i,2] - x[i,4]) ^ 4 for i in 1:length(x[:,1])]
    return f
end

function alphas(x)
    f = 1 / (-2 * b0 * (-0.3505067 - 0.5 * log(8 / x)))
    return f
end

function cutoff_true(x)
    f = [1 + x[i,1] / 8 * Lambda ^ 2 * 8 * pi * b0 * alphas(x[i,1]) ^ G1 + 0 * x[i,1] / 8 * Lambda ^ 2 * 8 * pi * b0 * alphas(x[i,1]) ^ G2 - 0 * (x[i,1] / 8) ^ 2 for i in 1:length(x[:,1])]
    return f
end

function model_true(x,p)
    return SU3_continuum(x,p) .* cutoff_true(x)
end

function ChPT_a2(x,p) 
    f = SU3_continuum(x,p) .* [1 + p[3] * x[i,1] / 8 for i in 1:length(x[:,1])]
    return f
end

function ChPT_a2alphas(x,p) 
    f = SU3_continuum(x,p) .* [1 + p[3] * x[i,1] / 8 * alphas(x[i,1]) ^ p[4] for i in 1:length(x[:,1])]
    return f
end

function ChPT_a2phi2(x,p) 
    f = SU3_continuum(x,p) .* [1 + p[3] * x[i,1] / 8 + p[4] * x[i,1] / 8 * x[i,2] for i in 1:length(x[:,1])]
    return f
end

function Taylor2_a2(x,p) 
    f = Tay_continuum(x,p) .* [1 + p[3] * x[i,1] / 8 for i in 1:length(x[:,1])]
    return f
end

function Taylor2_a2phi2(x,p) 
    f = Tay_continuum(x,p) .* [1 + p[3] * x[i,1] / 8 + p[4] * x[i,1] / 8 * x[i,2] for i in 1:length(x[:,1])]
    return f
end

function Taylor4_a2(x,p) 
    f = Tay4_continuum(x,p) .* [1 + p[4] * x[i,1] / 8 for i in 1:length(x[:,1])]
    return f
end

function Taylor4_a4(x,p) 
    f = Tay4_continuum(x,p) .* [1 + p[4] * x[i,1] / 8 + p[5] * (x[i,1] / 8) ^ 2  for i in 1:length(x[:,1])]
    return f
end

################################################################

phi2beta1 = [0.078, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.74];
phi2beta2 = [0.08, 0.1, 0.15, 0.2, 0.3, 0.33, 0.4, 0.45, 0.5, 0.52, 0.74];
phi2beta3 = [0.078, 0.18, 0.25, 0.3, 0.35, 0.38, 0.4, 0.45, 0.52, 0.6, 0.65, 0.73];
phi2beta4 = [0.078, 0.1, 0.13, 0.17, 0.2, 0.25, 0.29, 0.51, 0.55, 0.6, 0.73];
phi2beta5 = [0.079, 0.1, 0.2, 0.22, 0.25, 0.3, 0.33, 0.35, 0.5, 0.55, 0.6, 0.73];
t0 = [2.89, 3.66, 5.17, 8.6, 14.0];
x = [[[1 / t0[1] for i in 1:length(phi2beta1)] phi2beta1 [1.1 for i in 1:length(phi2beta1)] [0.73 for i in 1:length(phi2beta1)]]; 
     [[1 / t0[2] for i in 1:length(phi2beta2)] phi2beta2 [1.1 for i in 1:length(phi2beta2)] [0.73 for i in 1:length(phi2beta2)]];
     [[1 / t0[3] for i in 1:length(phi2beta3)] phi2beta3 [1.1 for i in 1:length(phi2beta3)] [0.73 for i in 1:length(phi2beta3)]];
     [[1 / t0[4] for i in 1:length(phi2beta4)] phi2beta4 [1.1 for i in 1:length(phi2beta4)] [0.73 for i in 1:length(phi2beta4)]];
     [[1 / t0[5] for i in 1:length(phi2beta5)] phi2beta5 [1.1 for i in 1:length(phi2beta5)] [0.73 for i in 1:length(phi2beta5)]]];

ptrue = [4.77,-0.45];
xph = [0.0 0.078 1.1 0.73];

n = length(x[:,1]);
dat = model_true(x,ptrue) .+ model_true(x,ptrue) .* [uwreal(randn(10000), i) for i in 1:n] ## 1% errors
dat_ref = deepcopy(dat)
uwerr.(dat);

#=
fig = figure("large")
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 20
xlabel(L"$\phi_2$")
ylabel(L"$\sqrt{t_0}f_{\pi K}$")
errorbar(x[:,2], value.(dat), err.(dat), label="", fmt="s", capsize=10.0)
ax = gca()
#ax[:set_ylim]([0.08, 0.115])
#legend()
tight_layout()
#close("all")
=#

cut = [collect(1:n), [collect(1:11); collect(13:22); collect(24:33); collect(35:44); collect(46:56)], [collect(1:8); collect(13:20); collect(24:31); collect(35:42); collect(46:54)], collect(13:n), [collect(13:22); collect(24:33); collect(35:44); collect(46:56)]]

y = [dat[cut[i]] for i in 1:length(cut)];
xcut = [x[cut[i],:] for i in 1:length(cut)];

models = [ChPT_a2; ChPT_a2phi2; ChPT_a2alphas; Taylor2_a2; Taylor2_a2phi2; Taylor4_a2; Taylor4_a4]
params = [3,4,4,3,4,4,5]

TIC = Array{Float64,1}()
pval_vec = Array{Float64,1}()
dof_vec = Array{Float64,1}()
yph = Array{uwreal,1}()

for j in 1:length(cut)
    for i in 1:length(models)
        upar, chi_exp, chi2, pval, doff = fit_alg(models[i], xcut[j], y[j], params[i])
        push!(TIC, chi2 - 2 * chi_exp)
        push!(dof_vec, doff)
        push!(yph, models[i](xph, upar)[1])
        push!(pval_vec, pval)
        println(i, " ", j)
    end
end
W = exp.(-0.5 .* TIC) ./ sum(exp.(-0.5 .* TIC));
W1 = deepcopy(W)
pv1 = deepcopy(pval_vec)
yph1 = deepcopy(yph)
yph_av = sum(yph .* W);
syst2 = sum(yph .^ 2 .* W) - yph_av ^ 2;
dof1 = deepcopy(dof_vec)

################################################################

yph_av = yph_av + uwreal([0.0,value(sqrt(syst2))], "syst")
uwerr(yph_av)
av = yph_av
uwerr.(yph)
W = W1

fig = figure("1",figsize=(5,5))
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 12

subplot(411) # Create the 2nd axis of a 3x1 arrax of axes
#mods = collect(1:length(sqrt_t0_vec))
x = collect(1:length(yph))
ax = gca() # Get the handle of the current axis
ylabel(L"$\sqrt{t_0}f_{\pi K}$")   
v = [value(av) for i in 1:length(yph)]
e = [err(av) for i in 1:length(yph)]
fill_between(x, v-e, v+e, color="deepskyblue", alpha=0.75)
errorbar(x, value.(yph), err.(yph), fmt="x", color="black")
#ax[:set_ylim]([0.139, 0.1477])
#ax[:set_xlim]([0, length(yph)+1])
#plt.xticks(x, mods)
setp(ax.get_xticklabels(),visible=false)

subplot(412) # Create the 3rd axis of a 3x1 array of axes
ax = gca()
#ax.set_yscale("log") # Set the x axis to a logarithmic scale
bar(x, pv1, color="deepskyblue")
#ax[:set_ylim]([0.0, 0.1])
ylabel("pval")
#plt.xticks(x, mods)
#xticks(rotation=90)
#ax[:set_xlim]([0, length(mods)+1])

subplot(413) # Create the 3rd axis of a 3x1 array of axes
ax = gca()
#ax.set_yscale("log") # Set the x axis to a logarithmic scale
bar(x, W, color="deepskyblue")
#ax[:set_ylim]([0.0, 0.1])
ylabel("W")
#plt.xticks(x, mods)
#xticks(rotation=90)
#ax[:set_xlim]([0, length(mods)+1])

subplot(414) # Create the 3rd axis of a 3x1 array of axes
ax = gca()
#ax.set_yscale("log") # Set the x axis to a logarithmic scale
bar(x, dof1, color="deepskyblue")
#ax[:set_ylim]([0.0, 0.1])
ylabel("dof")
#plt.xticks(x, mods)
#xticks(rotation=90)
#ax[:set_xlim]([0, length(mods)+1])

tight_layout()
subplots_adjust(hspace=0.15) 
